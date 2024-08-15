pacman::p_load(stars, sf, terra, data.table, tidyverse, lubridate,
               mgcv,
               viridis, scico, cols4all, 
               patchwork)


# Functions ===============================================
fn_dogliotti2015 <- function(p_w){
  # Dogliotti et al 2015, Remote Sensing of Environment
  # params for 645 nm
  # truncates negative vals and reflectance > 15%
  A_t <- 228.1
  C <- 0.164
  p_w <- ifelse(p_w < 0, 0, p_w)
  p_w <- ifelse(p_w > 0.15, 0.15, p_w)
  out <- A_t * (p_w) / 
    (1 - ((p_w)/ C))
  return(out)
}

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
# Torrens river discharge ------------------------
flow <- fread("data/daily_discharge_torrens_river_seaview_road_bridge.csv")
names(flow) <- c("date","time","dv")
flow[,`:=`(date = dmy(date))]
flow[,`:=`(year=year(date),month=month(date))]
flow <- flow[order(date)][,`:=`(dv_max5d = frollapply(dv,5,max,na.rm=T,align='right'))]
flow <- flow[order(date)][,`:=`(dv_max15d = frollapply(dv,15,max,na.rm=T,align='right'))]
flow <- flow[order(date)][,`:=`(dv_max30d = frollapply(dv,30,max,na.rm=T,align='right'))]
flow <- flow[order(date)][,`:=`(dv_sum5d = frollapply(dv,5,sum,na.rm=T,align='right'))]
flow <- flow[order(date)][,`:=`(dv_sum15d = frollapply(dv,15,sum,na.rm=T,align='right'))]
flow <- flow[order(date)][,`:=`(dv_sum30d = frollapply(dv,30,sum,na.rm=T,align='right'))]



## Data: MODIS-AQUA ocean color ==============================================
# kd490
aq <- terra::rast("data/processed_RS/Aqua_MODIS_rrs645.nc")
names(aq) <- time(aq)

daq <- aq %>% as.data.table(xy=T) %>% 
  pivot_longer(-c("x","y")) %>% 
  set_names(c("x","y","time","red")) %>% 
  as.data.table()

daq[,`:=`(dog = fn_dogliotti2015(red))]

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
daq[,`:=`(nobs = sum(is.na(dog)==F)),by=.(x,y)]
daq[,`:=`(nobs_zone_time = sum(is.na(dog)==F)),by=.(dz,time)]
daq[,`:=`(dog_u = mean(dog,na.rm=T),
          dog_sd = mean(dog,na.rm=T)),by=.(x,y,month)]
daq[,`:=`(dog_anom = dog - dog_u), by=.(x,y,month)]
daq[,`:=`(product = "AquaOceanColor")]


daq2 <- daq[is.na(dog_anom)==F][,.(
  val = weighted.mean(dog_anom, nobs_zone_time),
  val_sd = sd(dog_anom)),
  by=.(year,month,dz)][
    ,`:=`(d_ym = ymd(paste(year,month,1)))
  ][]
daq2 <- daq2[order(dz,d_ym)][,`:=`(val12 = frollmean(val,n=12,align='right')), 
                             by = dz]
daq2[,`:=`(dyear = decimal_date(d_ym))]
daq2[,`:=`(product = "AquaOceanColor")]



## Correlations =========================================
p1 <- daq2 %>% select(dyear, val, dz) %>% 
  pivot_wider(names_from = dz, 
              values_from = val) %>% 
  select(-dyear) %>% 
  drop_na() %>% 
  cor() %>% 
  ggcorrplot::ggcorrplot(lab=T, hc.order=F)+
  labs(title = "Monthly Dogliotti_2015 anom.")+
  scale_y_discrete(limits=rev); p1
p2 <- daq2 %>% select(dyear, val12, dz) %>% 
  pivot_wider(names_from = dz, 
              values_from = val12) %>% 
  select(-dyear) %>% 
  drop_na() %>% 
  cor() %>% 
  ggcorrplot::ggcorrplot(lab=T, hc.order=F)+
  labs(title = "12-month Dogliotti_2015 anom.")+
  scale_y_discrete(limits=rev); p2
p_out <- (p1|p2)
ggsave(p_out,
       filename = "figures/figure_AquaL3SMI_dogliotti_timeseries_correlation-by-zone.png",
       width = 24,
       height = 10,
       units='cm',
       dpi =350, 
       device = grDevices::png)


# daq2 %>% select(dyear, val12, dz) %>% 
#   pivot_wider(names_from = dz, 
#               values_from = val12) %>% 
#   select(-dyear) %>% 
#   drop_na() %>% 
#   cor() %>% 
#   corrplot::corrplot(method='number', 
#                      title = "Dogliotti_2015 12-month anomaly correlation")


flow2 <- flow[,.(
         dv_mean30d = mean(dv,na.rm=T),
         dv_sum30d = mean(dv_sum30d,na.rm=T),
         dv_max30d = mean(dv_max30d,na.rm=T)),by=.(year,month)]
merge(daq2, flow2 %>% select(year,month,starts_with('dv')), by=c("year","month"),
      allow.cartesian=T) %>% 
  group_by(dz) %>% 
  drop_na() %>% 
  summarize(rho1 = cor(val, dv_mean30d),
            rho2 = cor(val,dv_sum30d),
            rho3 = cor(val, dv_max30d))
merge(daq2, flow %>% select(year,month,starts_with('dv')), by=c("year","month"),
      allow.cartesian=T) %>% 
  group_by(dz) %>% 
  drop_na() %>% 
  summarize(rho1 = cor(val, dv),
            rho2 = cor(val,dv_sum30d),
            rho3 = cor(val, dv_max30d))


flow[dv_max30d >= max(dv_max30d,na.rm=T)]$date %>% min

events <- data.table(time = c(ymd("2005-11-08"),
                              ymd("2006-10-01"),
                              ymd("2010-09-05"),
                              ymd("2014-07-01"),
                              ymd("2019-08-10"), 
                              ymd("2016-10-01"), 
                              ymd("2022-12-01")),
                     event = c("2005 Flood",
                               "2006 Dredging",
                               "2010 End of Millenium drought",
                               "2014 Penrice Closure",
                               "2019 dredge", 
                               "2016 storm, high flow",
                               "2022/2023 Murray plume")) %>% 
  mutate(dyear=decimal_date(time))

daq2 %>% 
  ggplot(aes(dyear, val12))+
  geom_hline(aes(yintercept = 0),color='black') +
  geom_vline(data=events, 
             aes(xintercept=dyear), 
             color = 'black',
             lwd = 1) + 
  geom_vline(data=events, 
             aes(xintercept=dyear, 
                 color = event)) + 
  geom_line(color='black',
            aes(group=product),
            lwd = 1.5) + 
  geom_line() + 
  scale_color_discrete_c4a_cat(palette = 'okabe') + 
  # scale_color_brewer(type='qual',palette = 2) +
  labs(x = NULL, 
       y = "12 mo mean Dogliotti_2015 anomaly (FNU)", 
       color = NULL) + 
  # scale_x_date(date_breaks = "2 years", 
  #              date_labels = "%Y", 
  #              expand = c(0,0)) + 
  facet_wrap(~dz,scales='fixed',ncol=2) + 
  fn_theme() + 
  theme(legend.position = 'bottom', 
        legend.text = element_text(size = 10))


ggsave(filename=
         paste0("figures/figure_AquaL3SMI_dogliotti-smooothed-anom_timeseries_",Sys.Date(),".png"), 
       width = 30,
       height = 30,
       scale = 0.75,
       units = 'cm',
       dpi = 350)
