pacman::p_load(signal, stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork)
getDTthreads()

# Functions ================================
## Theme ============
fn_theme <- function (base_size = 12, base_family = "", base_line_size = base_size/22, 
                      base_rect_size = base_size/22){
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


# Main ========================================

## Torrens river discharge ------------------------
flow <- fread("data/daily_discharge_torrens_river_seaview_road_bridge.csv")
names(flow) <- c("date","time","dv")
flow[,`:=`(date = dmy(date))]
flow[,`:=`(year=year(date),month=month(date),doy = yday(date))]
flow <- flow[order(date)]
flow <- flow[order(date)][,`:=`(dv_max5d = frollapply(dv,5,max,na.rm=T,align='right'))]
flow <- flow[order(date)][,`:=`(dv_max15d = frollapply(dv,15,max,na.rm=T,align='right'))]
flow <- flow[order(date)][,`:=`(dv_max30d = frollapply(dv,30,max,na.rm=T,align='right'))]
flow <- flow[order(date)][,`:=`(dv_sum5d = frollapply(dv,5,sum,na.rm=T,align='right'))]
flow <- flow[order(date)][,`:=`(dv_sum15d = frollapply(dv,15,sum,na.rm=T,align='right'))]
flow <- flow[order(date)][,`:=`(dv_sum30d = frollapply(dv,30,sum,na.rm=T,align='right'))]

tmp <- flow[1:500,]
tmp <- tmp %>% mutate(dv_sum5d_sg = sgolay(dv_sum5d, n= 11))



tmp %>% 
  mutate(dv_sum5d_sg = sgolayfilt(.$dv_sum5d, p=5, n= 15)) %>% 
  mutate(dv_sum5d_diff = c(NA,diff(.$dv_sum5d))) %>% 
  ggplot(aes(date, dv))+
  geom_line()+
  geom_line(aes(date,dv_sum5d/5), col='blue')+
  geom_line(aes(date,dv_sum5d_diff/5), col='red')

# phenofit::PhenoDeriv(tmp$dv_sum5d, t=tmp$doy, der1=0,IsPlot=F,ana)
# 
# tmp %>% mutate(test = phenofit::PhenoDeriv(dv_sum5d, t=doy))


tmp[1:365,]$dv_sum5d %>% plot(type='l')
sgolayfilt(flow[1:365,]$dv_sum5d, n = 5) %>% lines(col='red')







pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork)
getDTthreads()

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
myd <- rast("data/processed_RS/Aqua_MODIS_kd490.nc")
names(myd) <- time(myd)

dat0 <- myd %>% 
  as.data.table(xy=T)

dat0 <- dat0 %>% 
  pivot_longer(-c('x','y')) %>% 
  set_names(c("x","y","time","kd490")) %>% 
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
dat <- dat[is.na(kd490)==F]

dat <- dat[,`:=`(ref_p95 = quantile(kd490,0.95), 
          ref_p50 = median(kd490)),
    by=.(x,y,month)]

dat[time<ymd('2015-01-01')] %>%
  ggplot(aes(time, 100*(kd490-ref_p50)/ref_p50,color=dz))+
  geom_hline(aes(yintercept = 0), col='grey') + 
  geom_smooth(method='gam', 
              formula = y~s(x,bs='ad',k=50)) + 
  scale_color_brewer(type='qual',palette = 2) + 
  facet_wrap(~dz, ncol=1, scales='free')+
  fn_theme()

dat[time>=ymd('2015-01-01')][time <= ymd("2024-01-01")] %>%
  ggplot(aes(time, 100*(kd490-ref_p50)/ref_p50,color=dz))+
  geom_hline(aes(yintercept = 0), col='grey') + 
  geom_smooth(method='gam', 
              formula = y~s(x,bs='ad',k=50)) + 
  scale_color_brewer(type='qual',palette = 2) + 
  facet_wrap(~dz, ncol=1, scales='free')+
  fn_theme()

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



dat2[time>=ymd('2019-01-01')][time <= ymd("2020-01-01")] %>%
  ggplot(aes(time, kd490,color=dz))+
  geom_smooth()
  geom_smooth(method='gam',
              formula = y~s(x),
              method.args=list(discrete=T)) 




dat[time>=ymd('2019-06-01')][time <= ymd("2020-01-01")]


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
  facet_wrap(~dz,scales='free',ncol=2)+
  fn_theme()+
  theme(legend.position = 'bottom')
ggsave(filename = 'figures/AQUA_MODIS_kd490_exponentialDecayFollowingDredge.png',
       width = 15,
       height = 15,
       units='cm',
       dpi=350,
       device = grDevices::png)



# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # AMC --------------------------------
# amc <- vect("outputs/AMC_10m_split.shp")
# zones <- rasterize(project(amc, tmp),tmp,field = "dz") %>% 
#   as.data.table(xy=T)
# 
# zones$dz %>% table %>% 
#   as.data.table() %>% 
#   knitr::kable()
# 
# dat <- merge(dat0, zones,by=c('x','y'))
# dat[,`:=`(dz = factor(dz,
#                       ordered = T,
#                       levels = c("NorthZone_gt10", "NorthZone_lt10", 
#                                  "CentralZone_gt10", "CentralZone_lt10",  
#                                  "SouthZone_gt10", "SouthZone_lt10"),
#                       labels = c("NorthZone_gt10", "NorthZone_lt10", 
#                                  "CentralZone_gt10", "CentralZone_lt10",  
#                                  "SouthZone_gt10", "SouthZone_lt10")))]
# dat[,`:=`(time = as_date(time))]
# dat[,`:=`(year=year(time),
#           month=month(time))]
# dat <- dat[is.na(kd490)==F]
# # dat <- dat %>% 
# #   mutate(kd490 = ifelse(kd490 > 2.5, 2.5, kd490))
# 
# dat[time>=ymd('2019-07-01')][time <= ymd("2020-01-01")] %>%
#   ggplot(aes(time, kd490,color=dz))+
#   geom_smooth(method='bam',
#               formula = y~s(x),
#               method.args=list(discrete=T))
# 
# m0 <- dat[time>=ymd('2019-07-01')][time <= ymd("2020-01-01")] %>% 
#   mutate(tsd = time-ymd("2019-07-01")) %>% 
#   mutate(tsd = as.numeric(tsd)) %>% 
#   # pull(tsd) %>% class
#   bam(kd490 ~ log(tsd), 
#       data=.) 
# 
# peak_time <- dat[time>=ymd('2019-06-01')][time <= ymd("2020-01-01")] %>%
#   .[,.(hi = quantile(kd490,0.9)),by=.(dz,time)] %>% 
#   group_by(dz) %>% 
#   summarize(event = min(time[hi==max(hi)]), 
#             val = max(hi)) %>% 
#   ungroup() %>% 
#   arrange(dz) %>% 
#   as.data.table()
# 
# tmp2 <- merge(peak_time, 
#               dat[time>=ymd('2019-06-01')][time <= ymd("2020-01-01")], 
#               by='dz')
# tmp2 <- tmp2 %>% filter(time >= event) %>% 
#   mutate(tsd = as.numeric(time - event))
# tmp2 %>% dim
# 
# vec_zones <- unique(tmp2$dz)
# 
# gam(kd490 ~ log(I(tsd+0.5)), 
#     data=tmp2[is.na(kd490)==F][dz==vec_zones[6]]) %>% 
#   plot(all.terms=T)
# 
# m3 <- bam(kd490 ~ te(x,y) + log(I(tsd+0.5))*dz,
#           data=tmp2,
#           discrete=T)
# summary(m3)
# tmp_pred <- expand_grid(unique(tmp2[,.(x,y,dz)]), tsd=c(1:30,seq(from=35,to=180,length.out=10))) %>% 
#   mutate(pred = predict(m3, newdata=.,type='response')) %>% 
#   as.data.table()
# 
# tmp_pred[,.(val = median(pred)),by=.(dz,tsd)] %>% 
#   filter(tsd <= 100) %>% 
#   ggplot(aes(tsd,val,color=dz))+
#   geom_line()+
#   scale_color_brewer(type='qual',palette = 2) + 
#   labs(x="Days since max kd490 during dredge",
#        y="Zonal median kd490",
#        title = 'Exponential decay following dredge event') +
#   # facet_wrap(~dz,scales='free')+
#   fn_theme()+
#   theme(legend.position = 'bottom')
# ggsave(filename = 'figures/S2_hack-kd490_exponentialDecayFollowingDredge.png',
#        width = 15,
#        height = 15,
#        units='cm',
#        dpi=350,
#        device = grDevices::png)
