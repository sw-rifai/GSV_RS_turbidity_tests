

# Functions ===============================================
fn_kd490_mod <- function(bg){
  # https://oceancolor.gsfc.nasa.gov/resources/atbd/kd/
  a0 <- -0.8813
  a1 <- -2.0584
  a2 <- 2.5878
  a3 <- -3.4885
  a4 <- -1.5061
  log10_kd490 <- a0 + a1*log10(bg) + a2*(log10(bg))**2 + a3*(log10(bg))**3 + a4*(log10(bg))**4
  kd490 <- 0.0166 + 10**log10_kd490
  kd490 <- ifelse(bg < 0.15, NA, kd490)
  return(kd490)
}

fn_kd490_oli <- function(bg){
  # https://oceancolor.gsfc.nasa.gov/resources/atbd/kd/
  vec_params <- c(-0.9054,	-1.5245,	2.2392,	-2.4777,	-1.1099)
  
  a0 <- vec_params[1]
  a1 <- vec_params[2]
  a2 <- vec_params[3]
  a3 <- vec_params[4]
  a4 <- vec_params[5]
  log10_kd490 <- a0 + a1*log10(bg) + a2*(log10(bg))**2 + a3*(log10(bg))**3 + a4*(log10(bg))**4
  kd490 <- 0.0166 + 10**log10_kd490
  return(kd490)
}

fn_theme <- function (base_size = 14, base_family = "", base_line_size = base_size/22, 
                      base_rect_size = base_size/22) {
  half_line <- base_size/2
  theme_bw(base_size = base_size, base_family = base_family, 
           base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace% 
    theme(axis.text = element_text(colour = "black", size = rel(0.8)), 
          axis.ticks = element_line(colour = "black", 
                                    linewidth = rel(0.5)),
          axis.ticks.length = unit(-0.1,'cm'),
          panel.border = element_rect(fill = NA, colour = "black", 
                                      linewidth = rel(1)), panel.grid = element_line(colour = "black"), 
          panel.grid.major = element_blank(), #element_line(linewidth = rel(0.1)), 
          panel.grid.minor = element_blank(), #element_line(linewidth = rel(0.05)), 
          strip.background = element_rect(fill = "transparent",
                                          color = 'transparent'), 
          strip.text = element_text(colour = "black", size = rel(0.8), 
                                    margin = margin(0.8 * half_line, 0.8 * half_line, 
                                                    0.8 * half_line, 0.8 * half_line)), 
          legend.position = c(0.5,0.025),
          legend.justification = c(0.5,0.025),
          legend.text = element_text(size = rel(0.8)),
          legend.direction = 'horizontal',
          legend.background =  element_rect(fill='transparent', 
                                            color='transparent'),
          complete = TRUE 
          # axis.text.y = element_blank()
          )
}




# Main ========================

data.table(bg = seq(0.15,0.25,length.out=1000)) %>% 
  mutate(MODIS = fn_kd490_mod(bg),
         L8_OLI = fn_kd490_oli(bg)) %>% 
  pivot_longer(-bg) %>% 
  ggplot(aes(bg,value,color=name))+
  geom_line(lwd=1)+
  scale_y_log10() +
  labs(x='Blue/Green',y='kd490',color="SeaDAS parameter set") + 
  cols4all::scale_color_discrete_c4a_cat(palette = "miscs.okabe") + 
  theme_linedraw()+
  theme(legend.position = 'bottom')


dat_m <- arrow::read_parquet("outputs/MCD43A4_AMC_blue_to_green.parquet") %>% 
  mutate(year = year(time),
         month = month(time), 
         doy = yday(time)) %>% 
  mutate(ddate = decimal_date(time)) %>% 
  mutate(fzone = factor(dz)) %>% 
  mutate(kd490_mean = fn_kd490_mod(mean))

# table(is.na(dat_m$kd490_mean))
quantile(dat_m$kd490_mean, c(0.005, 0.995),na.rm=T)
# dat_m$mean %>% hist(100); abline(v=0.15)
# dat_m$kd490_mean %>% hist(100)

# 'mean' Weighted MCD43 anomaly =========================
m_ref_mcd <- bam(kd490_mean ~ 
                   fzone + 
                   s(doy,fzone,bs=c("fs","cc"),k=6),
                 data=dat_m[kd490_mean > 0.01][kd490_mean < 6],
                 weights = count,
                 family = Gamma(link='log'),
                 select = T, 
                 discrete = T)
summary(m_ref_mcd)
draw(m_ref_mcd)

dat_m <- dat_m %>% 
  mutate(kd490_u = predict(m_ref_mcd,type='response',newdata=.)) %>% 
  mutate(kd490_mean_anom = kd490_mean - kd490_u)

### Do we have the correct number of observation?
nrow(dat_m) == length(unique(dat_m$time))*length(unique(dat_m$dz))

dat_m2 <- dat_m[is.na(kd490_mean)==F][,.(
  val = weighted.mean(kd490_mean_anom, count),
  val_sd = sd(kd490_mean_anom)),
  by=.(year,month,fzone)][
    ,`:=`(d_ym = ymd(paste(year,month,1)))
  ][]
dat_m2 <- dat_m2[order(fzone,d_ym)][,`:=`(val12 = frollmean(val,n=12,align='right')), 
                                    by = fzone]

dat_m2 %>% 
  filter(d_ym >= ymd("2002-12-01")) %>% 
  ggplot(aes(d_ym, val12,color=fzone))+
  geom_hline(aes(yintercept = 0),color='black') +
  geom_line(color='black',
            aes(group=fzone),
            lwd = 1.5) + 
  geom_line() + 
  scale_color_brewer(type='qual',palette = 2) + 
  labs(x = NULL, 
       y = "12 mo mean kd490", 
       color = NULL) + 
  scale_x_date(date_breaks = "2 years", 
               date_labels = "%Y", 
               expand = c(0,0)) + 
  fn_theme() + 
  theme(legend.position = 'bottom')
ggsave(filename = paste0('figures/mcd43a4_weighted_kd490-zonalMean-anom_12mo_timeseries_',Sys.Date(),'.png'),
       width = 30,
       height = 15,
       scale = 0.75,
       units = 'cm',
       dpi = 350)
