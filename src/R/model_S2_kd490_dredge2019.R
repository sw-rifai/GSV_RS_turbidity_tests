pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork)
getDTthreads()

## Theme ============
fn_theme <- function (base_size = 11, base_family = "", base_line_size = base_size/22, 
                      base_rect_size = base_size/22) 
{
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

oz_poly <- rnaturalearth::ne_countries(country = "Australia",
                                       returnclass = "sf",
                                       scale = 'large')

# Torrens river discharge ------------------------
flow <- fread("data/daily_discharge_torrens_river_seaview_road_bridge.csv")
names(flow) <- c("date","time","dv")
flow[,`:=`(date = dmy(date))]
flow[,`:=`(year=year(date),month=month(date))]
flow <- flow[order(date)][,`:=`(dv_max5d = frollapply(dv,5,max,na.rm=T,align='right'))]
flow <- flow[order(date)][,`:=`(dv_max15d = frollapply(dv,15,max,na.rm=T,align='right'))]
flow <- flow[order(date)][,`:=`(dv_max30d = frollapply(dv,30,max,na.rm=T,align='right'))]


# kd490
s2 <- rast("data/processed_RS/S2-composited_hack-kd490_AMC_2019-01-02_2024-03-11.nc")
names(s2) <- time(s2)

tmp <- aggregate(s2,3,mean,na.rm=T,cores=12)
names(tmp)
time(tmp) <- time(s2)

dat0 <- tmp %>% as.data.table(xy=T) %>% 
  pivot_longer(-c('x','y')) %>% 
  set_names(c("x","y","time","kd490")) %>% 
  as.data.table()
dat0 <- dat0[is.na(kd490)==F]

dat0[,`:=`(time = as_date(time))]

# dat0[time<ymd("2024-01-10")][sample(.N,1e6)] %>% 
#   ggplot(aes(time, kd490))+
#   geom_smooth(method='bam',
#               formula = y~s(x),
#               method.args=list(discrete=T))
# dat0[time==ymd("2019-08-10")] %>% 
#   ggplot(aes(x,y,fill=kd490))+
#   geom_raster()+scale_fill_viridis_c()+coord_sf()

# AMC --------------------------------
amc <- vect("outputs/AMC_10m_split.shp")
zones <- rasterize(project(amc, tmp),tmp,field = "dz") %>% 
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
dat <- dat[is.na(kd490)==F]
dat <- dat %>% 
  mutate(kd490 = ifelse(kd490 > 2.5, 2.5, kd490))

dat[time>=ymd('2019-07-01')][time <= ymd("2020-01-01")] %>%
  ggplot(aes(time, kd490,color=dz))+
  geom_smooth(method='bam',
              formula = y~s(x),
              method.args=list(discrete=T))

m0 <- dat[time>=ymd('2019-07-01')][time <= ymd("2020-01-01")] %>% 
  mutate(tsd = time-ymd("2019-07-01")) %>% 
  mutate(tsd = as.numeric(tsd)) %>% 
  # pull(tsd) %>% class
  bam(kd490 ~ log(tsd), 
      data=.) 

peak_time <- dat[time>=ymd('2019-06-01')][time <= ymd("2020-01-01")] %>%
  .[,.(hi = quantile(kd490,0.9)),by=.(dz,time)] %>% 
  group_by(dz) %>% 
  summarize(event = min(time[hi==max(hi)]), 
            val = max(hi)) %>% 
  ungroup() %>% 
  arrange(dz) %>% 
  as.data.table()

tmp2 <- merge(peak_time, 
      dat[time>=ymd('2019-06-01')][time <= ymd("2020-01-01")], 
      by='dz')
tmp2 <- tmp2 %>% filter(time >= event) %>% 
  mutate(tsd = as.numeric(time - event))
tmp2 %>% dim

vec_zones <- unique(tmp2$dz)

gam(kd490 ~ log(I(tsd+0.5)), 
    data=tmp2[is.na(kd490)==F][dz==vec_zones[6]]) %>% 
  plot(all.terms=T)

m3 <- bam(kd490 ~ te(x,y) + log(I(tsd+0.5))*dz,
           data=tmp2,
          discrete=T)
summary(m3)
tmp_pred <- expand_grid(unique(tmp2[,.(x,y,dz)]), tsd=c(1:30,seq(from=35,to=180,length.out=10))) %>% 
  mutate(pred = predict(m3, newdata=.,type='response')) %>% 
  as.data.table()

tmp_pred[,.(val = median(pred)),by=.(dz,tsd)] %>% 
  filter(tsd <= 100) %>% 
  ggplot(aes(tsd,val,color=dz))+
  geom_line()+
  scale_color_brewer(type='qual',palette = 2) + 
  labs(x="Days since max kd490 during dredge",
       y="Zonal median kd490",
       title = 'Exponential decay following dredge event') +
  # facet_wrap(~dz,scales='free')+
  fn_theme()+
  theme(legend.position = 'bottom')
ggsave(filename = 'figures/S2_hack-kd490_exponentialDecayFollowingDredge.png',
       width = 15,
       height = 15,
       units='cm',
       dpi=350,
       device = grDevices::png)

getViz(m3) %>% 
  plot(allTerms=T) %>% 
  print(pages=1)

fit_decay <- function(zone){
  out <- gam(kd490 ~ te(x,y) + log(I(tsd+0.5)),
             data=tmp2 %>% filter(dz==zone))
  out$dz <- zone
  return(out)
}

fits <- lapply(vec_zones, fit_decay)

fits %>% lapply(., summary)


fits %>% lapply(., 
  FUN = function(m){
  pdat <- data.table(tsd = 1:150) %>% 
    mutate(pred = predict(m, newdata=.,type='response'),
           dz = m$dz)
  return(pdat)
}) %>% 
  rbindlist(.) %>% 
  ggplot(aes(tsd, pred,color=dz))+
  geom_line(aes(group=dz),color='black',lwd=1.25) + 
  geom_line() + 
  scale_color_brewer(type='qual',palette = 2) + 
  facet_wrap(~dz,scales='free')



library(nls.multstart)
fit_nls_decay_p3 <- function(zone){
  sm <- tmp2[dz==zone][,.(kd490 = quantile(kd490,0.99)),by=tsd]
    
  out <- nls_multstart(kd490 ~ b0 + b1*log(b2*(tsd+0.5)),
             data=sm,
             iter = 1000,
             supp_errors = 'Y',
             start_lower = c('b0'=-1,
                             'b1'=-1,
                             'b2'=-1),
             start_upper = c('b0'=1,
                             'b1'=1,
                             'b2'=1))
  # out$dz <- zone
  return(out)
}
fit_nls_decay_p2 <- function(zone){
  sm <- tmp2[dz==zone][,.(kd490 = quantile(kd490,0.99)),by=tsd]
  
  out <- nls_multstart(kd490 ~ b0 + b2*log((tsd+0.5)),
                       data=sm,
                       iter = 1000,
                       supp_errors = 'Y',
                       start_lower = c('b0'=-1,
                                       # 'b1'=-1,
                                       'b2'=-1),
                       start_upper = c('b0'=1,
                                       # 'b1'=1,
                                       'b2'=1))
  # out$dz <- zone
  return(out)
}
m2 <- fit_nls_decay_p2(vec_zones[1])
m3 <- fit_nls_decay_p3(vec_zones[1])
bbmle::AICtab(m2,m3,fits[[1]])



res_decay

fn <- function(x){
  predict(m0, newdata=data.table(tsd=x))
}
curve(fn(x),1,100)


ref <- dat[,.(k_p95 = quantile(kd490,0.95), 
              k_p50 = median(kd490)),by=.(x,y,month)]


# ref %>% 
#   ggplot(aes(x,y,fill=k_p95))+
#   geom_sf(data = oz_poly, 
#           inherit.aes = F) + 
#   geom_raster()+scale_fill_viridis_c(option='H')+
#   coord_sf(xlim = c(range(dat$x)),
#            ylim=range(dat$y),
#            expand = F)+
#   facet_wrap(~month,nrow=2) + 
#   labs(x=NULL,
#        y=NULL,
#        fill="kd490 \n95pct") + 
#   fn_theme() + 
#   theme(axis.text.x = element_blank(), 
#         legend.position = 'right',
#         legend.direction = 'vertical')
# ggsave(filename='figures/Aqua-kd490_ref-p95-monthly.png',
#        width = 15,
#        height = 15,
#        units='cm',
#        dpi = 350,
#        scale=1.5)


dat2  <- merge(dat,ref,by=c('x','y','month'))

dat2[is.na(kd490)==F][,.(val = mean(kd490 - k_p50)),
                      by=.(dz,time)] %>% 
  mutate(year = year(time)) %>% 
  .[,.(val = max(val)),by=.(year,dz)] %>% 
  filter(year < 2024) %>% 
  ggplot(aes(year, val))+
  geom_point() + 
  geom_smooth(method='lm',se=F) + 
  geom_line(lty=2)+
  labs(y = "Annual max kd490 anom.",
       x = NULL) + 
  facet_wrap(~dz,nrow=3) +
  fn_theme()
ggsave(filename='figures/Aqua-kd490_AMC-max_annual_anom.png',
       width = 15,
       height = 20,
       units='cm',
       dpi = 350,
       scale=1)


dat2[is.na(kd490)==F][,.(val = mean(kd490 - k_p50)),
                      by=.(dz,time)] %>% 
  mutate(year = year(time)) %>% 
  .[,.(val = mean(val)),by=.(year,dz)] %>% 
  filter(year < 2024) %>% 
  ggplot(aes(year, val))+
  geom_point() + 
  geom_smooth(method='lm',se=F) + 
  geom_line(lty=2)+
  labs(y = "Annual mean kd490 anom.",
       x = NULL) + 
  facet_wrap(~dz,nrow=3) +
  fn_theme()
ggsave(filename='figures/Aqua-kd490_AMC-mean_annual_anom.png',
       width = 15,
       height = 20,
       units='cm',
       dpi = 350,
       scale=1)


dat2[is.na(kd490)==F][,.(val = mean(kd490 - k_p50)),
                      by=.(dz,time)] %>% 
  mutate(year = year(time)) %>% 
  .[,.(val = min(val)),by=.(year,dz)] %>% 
  filter(year < 2024) %>% 
  ggplot(aes(year, val))+
  geom_point() + 
  geom_smooth(method='lm',se=F) + 
  geom_line(lty=2)+
  labs(y = "Annual min kd490 anom.",
       x = NULL) + 
  facet_wrap(~dz,nrow=3) +
  fn_theme()
ggsave(filename='figures/Aqua-kd490_AMC-min_annual_anom.png',
       width = 15,
       height = 20,
       units='cm',
       dpi = 350,
       scale=1)


dat2[,.(val = max(kd490,na.rm=T)),by=.(time,dz)] %>% 
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
ggsave(filename='figures/Aqua-kd490_anomaliesByZone.png',
       width = 25,
       height = 20,
       units='cm',
       device = grDevices::png,
       dpi = 350,
       scale=1)


dat2[is.na(kd490)==F][,.(val = mean(kd490 - k_p50)),
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

b1 <- bam(kd490 ~ 
            s(month,bs='cc') + 
            s(x,y,by=ddate),
          family = Gamma(link='log'),
          discrete = T,
          select = T,
          data=dat0[is.na(kd490)==F][time <= max_time][kd490>0] %>% 
            mutate(log10_kd490 = log10(kd490)))
summary(b1)
hist(dat0$kd490,100)

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
       title = "MODIS Aqua kd490 annual trend",
       fill=expression(paste(frac(paste(Delta,"kd490"),paste(Delta,'year'))))
  )+ 
  coord_sf(xlim = c(range(dat$x)),
           ylim=range(dat$y),
           expand = F) + 
  fn_theme() +
  theme(axis.text.x = element_blank(), 
        legend.position = 'right',
        legend.direction = 'vertical',
        text = element_text(family = 'Helvetica'))
ggsave(filename='figures/Aqua-kd490_annual_trend.png',
       width = 15,
       height = 20,
       units='cm',
       device = grDevices::png,
       dpi = 350,
       scale=1)





## 

b2 <- bam(kd490 ~ 
            te(month,year,dz,bs=c('cc','tp','fs'),k=10), 
          family = Gamma(link='log'),
          discrete = T,
          select = T,
          data=dat2[is.na(kd490)==F])
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
ggsave(filename='figures/Aqua-kd490_seasonalCycleByZone.png',
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
#        fill=expression(paste(frac(paste(Delta,"kd490"),paste(Delta,'year'))))
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
# bam(kd490 ~ year, 
#     family = Gamma(link='log'),
#     data=dat0) %>% 
#   summary()
# 
