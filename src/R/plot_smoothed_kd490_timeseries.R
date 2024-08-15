pacman::p_load(stars, sf, terra, data.table, tidyverse, lubridate,
               mgcv,
               viridis, scico, cols4all, 
               patchwork)


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

# Main ==================================


## Data: MCD43A4 -----------------------
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

### mean Weighted MCD43 anomaly =========================
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
# ggsave(filename = paste0('figures/mcd43a4_weighted_kd490-zonalMean-anom_12mo_timeseries_',Sys.Date(),'.png'),
#        width = 30,
#        height = 15,
#        scale = 0.75,
#        units = 'cm',
#        dpi = 350)


## Data: MODIS-AQUA ocean color ==============================================
# kd490
aq <- terra::rast("data/processed_RS/Aqua_MODIS_kd490.nc")
names(aq) <- time(aq)

daq <- aq %>% as.data.table(xy=T) %>% 
  pivot_longer(-c("x","y")) %>% 
  set_names(c("x","y","time","kd")) %>% 
  as.data.table()

# AMC 
amc <- vect("outputs/AMC_10m_split.shp")
zones <- rasterize(project(amc, aq),aq,field = "dz") %>% 
  as.data.table(xy=T)

zones$dz %>% table %>% 
  as.data.table() %>% 
  knitr::kable()

daq <- merge(daq, zones,by=c('x','y'))
daq[,`:=`(dz = factor(dz,
                      ordered = T,
                      levels = c("NorthZone_gt10", "NorthZone_lt10", 
                                 "CentralZone_gt10", "CentralZone_lt10",  
                                 "SouthZone_gt10", "SouthZone_lt10"),
                      labels = c("NorthZone_gt10", "NorthZone_lt10", 
                                 "CentralZone_gt10", "CentralZone_lt10",  
                                 "SouthZone_gt10", "SouthZone_lt10")))]


daq[,`:=`(time = as_date(time))]
daq[,`:=`(dyear = decimal_date(time))]
daq[,`:=`(dyear_c = dyear - min(daq$dyear))]
daq[,`:=`(year=year(time),
         month=month(time))]
daq[,`:=`(nobs = sum(is.na(kd)==F)),by=.(x,y)]
daq[,`:=`(nobs_zone_time = sum(is.na(kd)==F)),by=.(dz,time)]
daq[,`:=`(kd_u = mean(kd,na.rm=T),
         kd_sd = mean(kd,na.rm=T)),by=.(x,y,month)]
daq[,`:=`(kd_anom = kd - kd_u), by=.(x,y,month)]
daq[,`:=`(product = "AquaOceanColor")]


daq2 <- daq[is.na(kd_anom)==F][,.(
  val = weighted.mean(kd_anom, nobs_zone_time),
  val_sd = sd(kd_anom)),
  by=.(year,month,dz)][
    ,`:=`(d_ym = ymd(paste(year,month,1)))
  ][]
daq2 <- daq2[order(dz,d_ym)][,`:=`(val12 = frollmean(val,n=12,align='right')), 
                                    by = dz]
daq2[,`:=`(dyear = decimal_date(d_ym))]
daq2[,`:=`(product = "AquaOceanColor")]



### Data: MCD43A4
dmcd <- arrow::read_parquet("outputs/MCD43A4_AMC_blue_to_green.parquet") %>% 
  mutate(year = year(time),
         month = month(time), 
         doy = yday(time)) %>% 
  mutate(dyear = decimal_date(time)) %>% 
  mutate(dz = factor(dz)) %>% 
  mutate(kd = fn_kd490_mod(mean), 
         product = "MCD43A4")
dmcd[,`:=`(dz = factor(dz,
          ordered = T,
          levels = c("NorthZone_gt10", "NorthZone_lt10", 
                     "CentralZone_gt10", "CentralZone_lt10",  
                     "SouthZone_gt10", "SouthZone_lt10"),
          labels = c("NorthZone_gt10", "NorthZone_lt10", 
                     "CentralZone_gt10", "CentralZone_lt10",  
                     "SouthZone_gt10", "SouthZone_lt10")))]

dmcd[,`:=`(kd_u = mean(kd,na.rm=T),
          kd_sd = mean(kd,na.rm=T)),by=.(dz,month)]
dmcd[,`:=`(kd_anom = kd - kd_u), by=.(dz,month)]

dmcd2 <- dmcd[is.na(kd_anom)==F][,.(
  val = weighted.mean(kd_anom, count),
  val_sd = sd(kd_anom)),
  by=.(year,month,dz)][
    ,`:=`(d_ym = ymd(paste(year,month,1)))
  ][]
dmcd2 <- dmcd2[order(dz,d_ym)][,`:=`(val12 = frollmean(val,n=12,align='right')), 
                                    by = dz]
dmcd2[,`:=`(dyear = decimal_date(d_ym))]
dmcd2[,`:=`(product = "MCD43A4")]


bind_rows(
daq2[d_ym >=ymd("2002-12-01")][d_ym <= ymd("2023-10-01")] %>% 
  select(dyear,val12,dz,product),
dmcd2[d_ym >=ymd("2002-12-01")][d_ym <= ymd("2023-10-01")] %>% 
  select(dyear,val12,dz,product)) %>% 
  ggplot(aes(dyear, val12,color=product))+
  geom_hline(aes(yintercept = 0),color='black') +
  geom_line(color='black',
            aes(group=product),
            lwd = 1.5) + 
  geom_line() + 
  scale_color_brewer(type='qual',palette = 2) +
  labs(x = NULL, 
       y = "12 mo mean kd490 anomaly", 
       color = NULL) + 
  # scale_x_date(date_breaks = "2 years", 
  #              date_labels = "%Y", 
  #              expand = c(0,0)) + 
  facet_wrap(~dz,scales='free') + 
  fn_theme() + 
  theme(legend.position = 'bottom')

ggsave(filename=
         paste0("figures/plot_smoothed_kd490_timeseries_",Sys.Date(),".png"), 
       width = 30,
       height = 15,
       scale = 0.75,
       units = 'cm',
       dpi = 350)


# daq2 %>% 
#   filter(d_ym >= ymd("2002-12-01")) %>% 
#   ggplot(aes(d_ym, val12,color=dz))+
#   geom_hline(aes(yintercept = 0),color='black') +
#   geom_line(color='black',
#             aes(group=dz),
#             lwd = 1.5) + 
#   geom_line() + 
#   scale_color_brewer(type='qual',palette = 2) +
#   labs(x = NULL, 
#        y = "12 mo mean kd490 anomaly", 
#        color = NULL) + 
#   scale_x_date(date_breaks = "2 years", 
#                date_labels = "%Y", 
#                expand = c(0,0)) + 
#   facet_wrap(~dz,scales='free') + 
#   fn_theme() + 
#   theme(legend.position = 'bottom')
# 
# 
# 
# bind_rows(
#   daq2[dz=="NorthZone_gt10"][dyear %between% c(2015,2019)] %>% 
#     select(dyear,val12,product),
#   dmcd2[dz=="NorthZone_gt10"][dyear %between% c(2015,2019)] %>% 
#     select(dyear,val12,product)) %>% 
#   ggplot(aes(dyear, val12,color=product))+
#   geom_line()
# 
# 
# 
# 
# 
# bind_rows(
#   daq[dz=="NorthZone_gt10"][dyear %between% c(2015,2019)] %>% 
#     select(dyear,kd,product),
#   dmcd[dz=="NorthZone_gt10"][dyear %between% c(2015,2019)] %>% 
#     select(dyear,kd,product)) %>% 
#   ggplot(aes(dyear, kd,color=product))+
#   geom_line()
