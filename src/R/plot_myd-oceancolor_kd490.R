pacman::p_load(terra,
               tidyverse, 
               data.table, 
               lubridate, 
               sf, 
               tidyterra, 
               viridis, 
               mgcv, 
               gratia)
# library(rasterVis)

amc <- vect("data/AMC/AnalysisExtentMask_GDA94z54.shp")
amc <- project(amc, crs("EPSG:4326"))
amc$Id <- c("central","north","south")
amc

myd <- terra::rast("data/processed_RS/Aqua_MODIS_kd490.nc")
terra:::readAll(myd)
terra::inMemory(myd)
time(myd)
# names(myd) <- time(myd)

tmp <- (myd[[time(myd) >= ymd("2022-09-01")]])
tmp <- tmp[[time(tmp) <= ymd("2023-10-24")]]
tmp <- log10(tmp)
names(tmp) <- time(tmp)

p_out <- ggplot()+geom_spatraster(data=tmp)+
  scale_fill_viridis_c(option="B")+
  coord_sf(expand = F) + 
  scale_x_continuous(breaks = c(138, 139)) + 
  facet_wrap(~lyr, nrow = 5) + 
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        legend.position = 'bottom',
        legend.key.width = unit(0.1,"npc"))
ggsave(p_out, 
       filename = "figures/myd_oceancolor_kd490_murray_flood.png", 
       width = 35, 
       height = 30, 
       units= 'cm',
       dpi = 350)



ragg::agg_png(filename = "figures/myd_oceancolor_kd490_murray_flood.png",
    width = 4000,
    height = 4000,
    units='px',
    pointsize = 36)
panel(log10(myd[[time(myd) >= ymd("2022-06-01")]][[time(myd) <= ymd("2023-12-31")]]), 
     col = inferno(30), 
     range = c(-1.2,0.6), 
     maxcell = 1e6,
     mar = c(0,0,0,0),
     # bltr
     # mar = c(1,0,2,0),
     axes = T,
     background = "grey70", 
     nc = 10,
     nr = 8,
     maxnl = 100,
     # par = list( mar = c(1,0,2,0)),
     breakby = "equint")
dev.off()

tmp <- log10(myd[[time(myd) >= ymd("2022-06-01")]][[time(myd) <= ymd("2023-12-31")]])
rasterVis::levelplot(tmp)


## 
dat <- myd %>% as.data.table(xy=T) %>% 
  pivot_longer(-c("x","y")) %>% 
  rename(time = name, 
         kd490 = value) %>% 
  mutate(date = lubridate::date(time))
dat <- dat %>% as.data.table()
dat[,.(x,y)] %>% unique
dat <- dat[,`:=`(year=year(time),
                 month = month(time),
                 doy = yday(time))]
dat <- dat[,`:=`(ddate = decimal_date(date))]

plot(mean(myd,na.rm=T))
lines(amc)

dat$kd490 %>% log10 %>% hist(1000)
m0 <- bam(log10(kd490) ~ s(x,y,bs='sos'), 
          data=dat, 
          discrete = T,
          select = T)
summary(m0)
plot(m0, scheme=4)


m1 <- bam(log10(kd490) ~ 
            s(doy,bs='cc') + 
            s(x,y,bs='sos'), 
          data=dat, 
          discrete = T,
          select = T)
summary(m1)
plot(m1,select=1)


m2 <- bam(log10(kd490) ~ 
            # s(ddate,k=20) + 
            te(x,y,doy, bs=c("tp","tp","cc")),# + 
            # s(x,y,bs='sos'), 
          data=dat, 
          discrete = T,
          select = T)
summary(m2)
plot(m2,select=1)

m2_s <- smooth_estimates(m2)

m2_s %>% 
  arrange(x) %>% 
  ggplot(aes(doy, .estimate, group=paste(x,y),color=y)) + 
  geom_line()+
  scale_color_viridis_c(option='H')


m3 <- bam(log10(kd490) ~ 
            s(ddate,k=20) +
            te(x,y,doy, bs=c("tp","tp","cc")),# + 
          # s(x,y,bs='sos'), 
          data=dat, 
          discrete = T,
          select = T)
summary(m3)
plot(m3,select=1)


m4 <- bam(log10(kd490) ~ 
            # te(x,y,ddate, k=c(10,10,20)) +
            s(doy,bs='cc') + 
            te(x,y,ddate, bs=c("tp","tp","tp")),# + 
          # s(x,y,bs='sos'), 
          data=dat, 
          discrete = T,
          select = T)
summary(m4)
plot(m4,select=1)
m4_s <- smooth_estimates(m4, n_3d=100) %>% setDT()

m4_s %>% select(x,y,doy, ddate) %>% unique
m4_s[is.na(ddate)==F][sample(.N, 1000)] %>% 
  ggplot(aes(ddate, .estimate))+
  geom_point()


m4_s[is.na(ddate)==F] %>% 
  mutate(xc = cut_interval(x,5), 
         yc= cut_interval(y, 5)) %>% 
  ggplot(aes(ddate, .estimate))+
  geom_smooth() +  
  facet_grid(factor(rev(yc), levels = rev(levels(yc))) ~ xc, scales = 'free')


## 
tmp <- terra::extract(myd, amc, xy=T, fun='mean', 
                      bind = T)

tmp2 <- tmp %>% 
  as.data.table() %>% 
  pivot_longer(-c("Id","gridcode","zone")) %>% 
  as.data.table()
tmp2 <- tmp2 %>%
  rename(time = name, 
         kd490 = value)
tmp2 <- tmp2 %>% 
  mutate(time = ymd(time))
tmp2 <- tmp2[,`:=`(year=year(time),
           month=month(time),
           doy=yday(time))]


tmp2[year==2019] %>% 
  ggplot(aes(time, kd490, color=zone))+
  geom_line()


tmp2 %>% 
  ggplot(aes(time, kd490, color=zone))+
  geom_line()


mean(myd,na.rm=T) %>% plot
terra::stdev(myd,na.rm=T) %>% plot
(terra::stdev(myd,na.rm=T)/mean(myd,na.rm=T)) %>% 
  log10() %>% 
  plot(col = inferno(100))



