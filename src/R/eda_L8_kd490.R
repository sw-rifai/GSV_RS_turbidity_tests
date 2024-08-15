pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork)
getDTthreads()

# Torrens river discharge ------------------------
flow <- fread("data/daily_discharge_torrens_river_seaview_road_bridge.csv")
names(flow) <- c("date","time","dv")
flow[,`:=`(date = dmy(date))]
flow[,`:=`(year=year(date),month=month(date))]
flow <- flow[order(date)][,`:=`(dv_max5d = frollapply(dv,5,max,na.rm=T,align='right'))]
flow <- flow[order(date)][,`:=`(dv_max15d = frollapply(dv,15,max,na.rm=T,align='right'))]
flow <- flow[order(date)][,`:=`(dv_max30d = frollapply(dv,30,max,na.rm=T,align='right'))]

## Load and collate raster stacks of Landsat
vec_dates <- fread("data/gee_SAWater_turbidity/DATES_LC08_C02_T1_L2_kd490_l8_2013-03-18_2023-09-18.csv")
ic <- rast("data/gee_SAWater_turbidity/LC08_C02_T1_L2_kd490_l8_2013-03-18_2023-09-18.tif")

# ## extract depth and AMC zones ------------------
amc <- vect("data/AMC/AnalysisExtentMask_GDA94z54.shp")
amc <- terra::project(amc, crs(ic))
amc <- amc %>% select(zone)


time(ic) <- vec_dates$date %>% 
  as_date()
ic <- crop(ic,amc)
ic_l8 <- clamp(ic,lower=0,upper=20)



val_u <- mean(ic_l8,na.rm=T)
val_max <- max(ic_l8,na.rm=T)
val_sd <- terra::stdev(ic_l8,na.rm=T)
val_count <- terra::countNA(ic_l8)
val_count <- (111 - val_count)
val_p95 <- terra::quantile(ic_l8, 0.95,na.rm=T)

plot(val_u)
plot(val_sd)
plot(val_count)
plot(val_max)
plot(val_p95)


tmp <- ic[[year(time(ic))==2019]]
names(tmp) <- time(tmp)
p1 <- ggplot()+ 
  geom_spatraster(data=tmp) + 
  geom_sf(data = amc, fill='transparent',
          color = 'white') + 
  scale_fill_viridis_c(option = 'B', 
                       # limits = c(-0.5,0.25),
                       limits = c(0,1),
                       oob =scales::squish) +
  coord_sf(expand = F) + 
  facet_wrap(~lyr, nrow = 2) + 
  theme_linedraw() + 
  theme(panel.grid=element_blank(), 
        axis.text.x = element_blank())
ggsave(plot = p1,
       filename = "figures/l8_kd490_2019.png",
       width = 20,
       height = 15,
       units = 'cm')


flow[year==2019] %>% 
  ggplot(aes(date,dv))+
  geom_line()
