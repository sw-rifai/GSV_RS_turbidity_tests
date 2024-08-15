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
vec_dates <- fread("data/gee_SAWater_turbidity/DATES_S2_kd490_s2_2019-09-01_2020-01-01.csv")
ss <- rast("data/gee_SAWater_turbidity/S2_kd490_s2_2019-09-01_2020-01-01.tif")
time(ss) <- vec_dates$date %>% as.Date()
names(ss) <- vec_dates$date %>% as.Date()
ss <- ss[[order(vec_dates$date)]]


# ## extract depth and AMC zones ------------------
amc <- vect("data/AMC/AnalysisExtentMask_GDA94z54.shp")
amc <- terra::project(amc, crs(ss))
amc <- amc %>% select(zone)


ss <- crop(ss,amc)
ss <- clamp(ss,lower=0,upper=20)



val_u <- mean(ss,na.rm=T)
val_max <- max(ss,na.rm=T)
val_sd <- terra::stdev(ss,na.rm=T)
val_count <- terra::countNA(ss)
val_count <- (111 - val_count)
val_p95 <- terra::quantile(ic_l8, 0.95,na.rm=T)

plot(val_u)
plot(val_sd)
plot(val_count)
plot(val_max)
plot(val_p95)



tmp2 <- ss[[year(time(ss))==2019]]
names(tmp2) <- time(tmp2)
p2 <- ggplot()+ 
  geom_spatraster(data=tmp2) + 
  geom_sf(data = amc, fill='transparent',
          color = 'white') + 
  scale_fill_viridis_c(option = 'B', 
                       # limits = c(-0.5,0.25),
                       limits = c(0.05,0.45),
                       oob =scales::squish) +
  coord_sf(expand = F) + 
  facet_wrap(~lyr, nrow = 2) + 
  theme_linedraw() + 
  theme(panel.grid=element_blank(), 
        axis.text.x = element_blank(), 
        strip.text = element_text(size=7))
ggsave(plot = p2,
       filename = "figures/S2_kd490_2019.png",
       width = 35,
       height = 15,
       units = 'cm')


# flow[year==2019] %>% 
#   ggplot(aes(date,dv))+
#   geom_line()
