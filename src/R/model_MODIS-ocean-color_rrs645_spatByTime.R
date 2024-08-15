pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork)
getDTthreads()


# Functions ---------------------------------------------------------------


## Theme ============
fn_theme <- function (base_size = 11, base_family = "", base_line_size = base_size/22, 
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
          legend.direction = 'horizontal',
          legend.background =  element_rect(fill='transparent', 
                                            color='transparent'),
          complete = TRUE)
}

fn_dogliotti <- function(p_w){
  # Dogliotti et al 2015, Remote Sensing of Environment
  # Notes: See ACOLITE for most recent parameter calibrations
  A_t <- 378.46
  C <- 0.19905
  out <- A_t * (p_w) / 
    (1 - ((p_w)/ C))
  return(out)
}




# Main --------------------------------------------------------------------

## Import data ====================

oz_poly <- rnaturalearth::ne_countries(country = "Australia",
                                       returnclass = "sf",
                                       scale = 'large')


# Torrens river discharge 
flow <- fread("data/daily_discharge_torrens_river_seaview_road_bridge.csv")
names(flow) <- c("date","time","dv")
flow[,`:=`(date = dmy(date))]
flow[,`:=`(year=year(date),month=month(date))]
flow <- flow[order(date)][,`:=`(dv_max5d = frollapply(dv,5,max,na.rm=T,align='right'))]
flow <- flow[order(date)][,`:=`(dv_max15d = frollapply(dv,15,max,na.rm=T,align='right'))]
flow <- flow[order(date)][,`:=`(dv_max30d = frollapply(dv,30,max,na.rm=T,align='right'))]


# AQUA nasa ocean color: dog
myd <- rast("data/processed_RS/Aqua_MODIS_rrs645.nc")
names(myd) <- time(myd)

## Process =====================
r_dog <- terra::app(myd, fn_dogliotti, cores = 12)



dat0 <- r_dog %>% 
  as.data.table(xy=T)

dat0 <- dat0 %>% 
  pivot_longer(-c('x','y')) %>% 
  set_names(c("x","y","time","dog")) %>% 
  as.data.table()
dat0[,`:=`(time = as_date(time))]

# AMC --------------------------------
amc <- vect("outputs/AMC_10m_split.shp")
zones <- rasterize(project(amc, myd),myd,field = "dz") %>% 
  as.data.table(xy=T)

zones$dz %>% table %>% 
  as.data.table() %>% 
  knitr::kable()

dat <- merge(dat0, zones,by=c('x','y'))
dat[,`:=`(dz = factor(dz,
                      ordered = T,
                      levels = c("NorthZone_gt10", "NorthZone_lt10", 
                                 "CentralZone_gt10", "CentralZone_lt10",  
                                 "SouthZone_gt10", "SouthZone_lt10"),
                      labels = c("NorthZone_gt10", "NorthZone_lt10", 
                                 "CentralZone_gt10", "CentralZone_lt10",  
                                 "SouthZone_gt10", "SouthZone_lt10")))]
dat[,`:=`(time = as_date(time))]
dat[,`:=`(year=year(time),
          month=month(time))]
dat <- dat[is.na(dog)==F]

dat <- dat[,`:=`(dog_p95 = quantile(dog,0.95), 
              dog_p50 = median(dog), 
              dog_u = mean(dog), 
              dog_sd = sd(dog)),
           by=.(x,y,month)]

dat2 <- dat

dat2[is.na(dog)==F][,.(val = mean(dog - dog_u)),
                      by=.(dz,time)] %>% 
  mutate(year = year(time)) %>% 
  .[,.(val = max(val)),by=.(year,dz)] %>% 
  filter(year < 2024) %>% 
  ggplot(aes(year, val))+
  geom_point() + 
  geom_smooth(method='lm',se=F) + 
  geom_line(lty=2)+
  labs(y = "Annual max dog anom.",
       x = NULL) + 
  facet_wrap(~dz,nrow=3) +
  fn_theme()
ggsave(filename='figures/Aqua-dog_AMC-max_annual_anom.png',
       width = 15,
       height = 20,
       units='cm',
       dpi = 350,
       scale=1)


dat2[is.na(dog)==F][,.(val = mean(dog - dog_u)),
                      by=.(dz,time)] %>% 
  mutate(year = year(time)) %>% 
  .[,.(val = mean(val)),by=.(year,dz)] %>% 
  filter(year < 2024) %>% 
  ggplot(aes(year, val))+
  geom_point() + 
  geom_smooth(method='lm',se=F) + 
  geom_line(lty=2)+
  labs(y = "Annual mean dog anom.",
       x = NULL) + 
  facet_wrap(~dz,nrow=3) +
  fn_theme()
ggsave(filename='figures/Aqua-dog_AMC-mean_annual_anom.png',
       width = 15,
       height = 20,
       units='cm',
       dpi = 350,
       scale=1)


dat2[is.na(dog)==F][,.(val = mean(dog - k_p50)),
                      by=.(dz,time)] %>% 
  mutate(year = year(time)) %>% 
  .[,.(val = min(val)),by=.(year,dz)] %>% 
  filter(year < 2024) %>% 
  ggplot(aes(year, val))+
  geom_point() + 
  geom_smooth(method='lm',se=F) + 
  geom_line(lty=2)+
  labs(y = "Annual min dog anom.",
       x = NULL) + 
  facet_wrap(~dz,nrow=3) +
  fn_theme()
ggsave(filename='figures/Aqua-dog_AMC-min_annual_anom.png',
       width = 15,
       height = 20,
       units='cm',
       dpi = 350,
       scale=1)


dat2[,.(val = max(dog,na.rm=T)),by=.(time,dz)] %>% 
  ggplot(aes(time, val))+
  geom_vline(aes(xintercept = ymd("2016-10-01")),
             lwd = 2,
             color='grey') +
  geom_vline(aes(xintercept = ymd("2019-08-10")),
             lwd = 2,
             color='grey') +
  geom_line(col="#cf0000")+
  geom_ribbon(data=dat2, 
              inherit.aes = F,
              aes(x=time,ymin=k_p50,ymax=k_p95),
              fill='grey30',
              alpha = 0.9,
              color='transparent',
              lty=0) + 
  facet_wrap(~dz, ncol=2,
             scales='free') + 
  fn_theme()
ggsave(filename='figures/Aqua-dog_anomaliesByZone.png',
       width = 25,
       height = 20,
       units='cm',
       device = grDevices::png,
       dpi = 350,
       scale=1)


dat2[is.na(dog)==F][,.(val = mean(dog - k_p50)),
                     by=.(dz,time)] %>% 
  ggplot(aes(time, val,color=dz))+
  geom_line()+
  facet_wrap(~dz,ncol=2,scales="free")



# GAM ===============================
library(gratia)
dat0[,`:=`(year=year(time),
           month = month(time))]
dat0[,`:=`(ddate = decimal_date(time))]
max_time <- (dat0$time %>% min) + years(22) + days(8)

b1 <- bam(dog ~ 
            s(month,bs='cc') + 
            s(x,y,by=ddate),
          family = Gamma(link='log'),
          discrete = T,
          select = T,
          data=dat0[is.na(dog)==F][time <= max_time][dog>0] %>% 
          mutate(log10_dog = log10(dog)))
summary(b1)
hist(dat0$dog,100)

residuals(b1) %>% hist(100)

expand_grid(unique(dat0[,.(x,y)]), 
            month = 7,
            ddate = c(2002.5,2023.5)) %>% 
  as.data.table() %>% 
  mutate(pred = predict(b1,newdata=., type='response')) %>% 
  group_by(x,y) %>% 
  summarize(delta = pred[ddate==2023.5] - pred[ddate==2002.5]) %>% 
  # pull(delta) %>% quantile(.,c(0.01,0.95))
  ungroup() %>% 
  ggplot(aes(x,y,fill=delta/(2023.5-2002.5)))+
  geom_raster()+
  geom_sf(data=oz_poly, 
          inherit.aes = F) + 
  scale_fill_viridis_c(option='B', 
                       limits = c(0,0.005),
                       oob=scales::squish) + 
  # scale_fill_continuous_c4a_div() + 
  labs(x=NULL,
       y=NULL,
       title = "MODIS Aqua dog annual trend",
       fill=expression(paste(frac(paste(Delta,"dog"),paste(Delta,'year'))))
  )+ 
  coord_sf(xlim = c(range(dat$x)),
           ylim=range(dat$y),
           expand = F) + 
  fn_theme() +
  theme(axis.text.x = element_blank(), 
        legend.position = 'right',
        legend.direction = 'vertical',
        text = element_text(family = 'Helvetica'))
ggsave(filename='figures/Aqua-dog_annual_trend.png',
       width = 15,
       height = 20,
       units='cm',
       device = grDevices::png,
       dpi = 350,
       scale=1)





## 

b2 <- bam(dog ~ 
            te(month,year,dz,bs=c('cc','tp','fs'),k=10), 
          family = Gamma(link='log'),
          discrete = T,
          select = T,
          data=dat2[is.na(dog)==F])
summary(b2)
tmp <- expand_grid(year=2002:2023,month=seq(1,to=12,length.out=100),
                   dz = factor(dat2$dz %>% levels)) %>% 
mutate(estimate =  predict(b2, newdata = .,
               type='response')) %>% 
  as.data.table()
tmp[,`:=`(dz = factor(dz,
                      ordered = T,
                      levels = c("NorthZone_gt10", "NorthZone_lt10", 
                                 "CentralZone_gt10", "CentralZone_lt10",  
                                 "SouthZone_gt10", "SouthZone_lt10"),
                      labels = c("NorthZone_gt10", "NorthZone_lt10", 
                                 "CentralZone_gt10", "CentralZone_lt10",  
                                 "SouthZone_gt10", "SouthZone_lt10")))]

tmp %>% 
  ggplot(aes(month, estimate,group=year,
             color=year))+
  geom_line() + 
  scale_color_viridis_c(option='H')+
  scale_x_continuous(breaks = 1:12) + 
  facet_wrap(~dz,ncol=2,scales='free_y') + 
  fn_theme()+
  theme( 
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.key.width = unit(0.1,'npc'),
        text = element_text(family = 'Helvetica'))
ggsave(filename='figures/Aqua-dog_seasonalCycleByZone.png',
       width = 15,
       height = 20,
       units='cm',
       device = grDevices::png,
       dpi = 350,
       scale=1)





# smooth_estimates(b1) %>% 
#   filter(is.na(x)==F) %>% 
#   ggplot(aes(x,y,fill=.estimate))+
#   geom_raster()+
#   geom_sf(data=oz_poly, 
#           inherit.aes = F) + 
#   scale_fill_viridis_c(option='B', 
#                        # limits = c(0.007,0.009),
#                        oob=scales::squish) + 
#   # scale_fill_continuous_c4a_div() + 
#   labs(x=NULL,
#        y=NULL,
#        fill=expression(paste(frac(paste(Delta,"dog"),paste(Delta,'year'))))
#        )+ 
#   coord_sf(xlim = c(range(dat$x)),
#            ylim=range(dat$y),
#            expand = F) + 
#   fn_theme() +
#   theme(axis.text.x = element_blank(), 
#         legend.position = 'right',
#         legend.direction = 'vertical')
# 
#   
# 
# bam(dog ~ year, 
#     family = Gamma(link='log'),
#     data=dat0) %>% 
#   summary()
# 
