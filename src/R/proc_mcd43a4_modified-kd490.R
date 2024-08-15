pacman::p_load(stars, sf, terra, data.table, tidyverse, lubridate,
               mgcv,
               viridis, scico, cols4all, 
               patchwork)

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

fn_kd490_mod <- function(bg){
  # https://oceancolor.gsfc.nasa.gov/resources/atbd/kd/
  
  # blue ~ 470
  # green ~ 555
  bg <- bg+0.2
  a0 <- -0.8813
  a1 <- -2.0584
  a2 <- 2.5878
  a3 <- -3.4885
  a4 <- -1.5061
  log10_kd490 <- a0 + a1*log10(bg) + a2*(log10(bg))**2 + a3*(log10(bg))**3 + a4*(log10(bg))**4
  kd490 <- 0.0166 + 10**log10_kd490
  kd490 <- ifelse(bg < 0.2, NA, kd490)
  return(kd490)
}
curve(fn_kd490_mod(x), 0.15, 2.96)

### DATA: MCD43A4
tmp <- terra::rast("data/processed_RS/MCD43A4-qaMasked_blue_to_green_AMC_monthly_2000-02-24_2023-10-22.nc")
terra:::readAll(tmp)
names(tmp) <- time(tmp)
tmp2 <- app(tmp, fun = function(i,ff) ff(i), 
            # cores=24,
            ff=fn_kd490_mod)
tmp2
time(tmp2) <- time(tmp)

tmp3 <- tmp2 %>% 
  as.data.table(xy=T) %>% 
  pivot_longer(-c("x","y")) %>% 
  set_names(c("x","y","time","kd")) %>% 
  as.data.table()


amc <- vect("outputs/AMC_10m_split.shp")
zones <- rasterize(project(amc, tmp),tmp,field = "dz") %>% 
  as.data.table(xy=T)

tmp4 <- merge(tmp3,zones,by=c("x","y"))

tmp4[,`:=`(dz = factor(dz,
                      ordered = T,
                      levels = c("NorthZone_gt10", "NorthZone_lt10", 
                                 "CentralZone_gt10", "CentralZone_lt10",  
                                 "SouthZone_gt10", "SouthZone_lt10"),
                      labels = c("NorthZone_gt10", "NorthZone_lt10", 
                                 "CentralZone_gt10", "CentralZone_lt10",  
                                 "SouthZone_gt10", "SouthZone_lt10")))]

dmcd <- tmp4 %>% 
  mutate(time = ymd(time)) %>% 
  mutate(dyear = decimal_date(time)) %>% 
  mutate(year=year(time),
         month=month(time))

dmcd <- dmcd[dyear >= 2002] # !!!

dmcd <- dmcd[,`:=`(kd_u = mean(kd,na.rm=T),
           kd_sd = mean(kd,na.rm=T)),
     by=.(x,y,dz,month)]
dmcd[,`:=`(kd_anom = kd - kd_u), by=.(x,y,dz,month)]
dmcd[,`:=`(count = sum(is.na(kd)==F)), 
     by=.(dz,dyear)]


dmcd2 <- dmcd[is.na(kd_anom)==F][,.(
  val = weighted.mean(kd_anom, count),
  val_sd = sd(kd_anom)),
  by=.(dyear,year,month,dz)][
    ,`:=`(d_ym = ymd(paste(year,month,1)))
  ]


dmcd2 <- dmcd2[order(dz,d_ym)][,`:=`(val12 = frollmean(val,n=12,align='right')), 
                               by = dz]
dmcd2[,`:=`(dyear = decimal_date(d_ym))]
dmcd2[,`:=`(product = "MCD43A4")]

dmcd2 %>% 
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
