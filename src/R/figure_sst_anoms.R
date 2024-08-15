pacman::p_load(terra,
               tidyverse, 
               data.table, 
               lubridate, 
               sf, 
               tidyterra) 

sst <- rast("data/processed_RS/Aqua_MODIS_sst_daytime.nc")
names(sst) <- time(sst)
amc <- vect("data/AMC/AnalysisExtentMask_GDA94z54.shp")
amc <- project(amc, crs("EPSG:4326"))

tmp <- terra::extract(sst, amc, fun='mean', xy=T)

dat <- tmp %>% 
  pivot_longer(-c("ID")) %>% 
  as.data.table() %>% 
  set_names(c("id", "date","sst"))

dat[,`:=`(time = ymd(date))]
dat[,`:=`(year=year(time),
          month=month(time),
          doy = yday(time))]
dat[,`:=`(sst_u = mean(sst,na.rm=T)), by=month]


dat[id==2] %>% 
  ggplot(aes(time, sst-sst_u,color=id))+
  geom_point()+
  geom_smooth(method='gam', 
              formula = y~s(x,bs='ad'))


dat[id==2][year %in% c(2003:2006, 2020:2024)][,.(anom = mean(sst-sst_u, na.rm=T)), by=.(id,year,month)] %>% 
 ggplot(aes(month, anom,color=factor(year),group=year))+
  geom_line(color='black',lwd=1.2) + 
  geom_line()+
  scale_color_viridis_d(option='H')+
  theme_linedraw()

dat[id==2][year >= 2020][,.(anom = mean(sst-sst_u, na.rm=T)), by=.(id,year,month)] 
  ggplot(aes(doy, sst-sst_u,color=year,group=year))+
  geom_line(color='black',lwd=1.2) + 
  geom_line()+
  scale_color_viridis_c(option='H')+
  theme_linedraw()
