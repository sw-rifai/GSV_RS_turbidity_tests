pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork)

scale_factor <- 0.0001;

# Import ============ 

# Torrens river discharge ------------------------
flow <- fread("data/daily_discharge_torrens_river_seaview_road_bridge.csv")
names(flow) <- c("date","time","dv")
flow[,`:=`(date = dmy(date))]
flow[,`:=`(year=year(date),month=month(date))]
flow <- flow[order(date)][,`:=`(dv_max5d = frollapply(dv,5,max,na.rm=T,align='right'))]
flow <- flow[order(date)][,`:=`(dv_max15d = frollapply(dv,15,max,na.rm=T,align='right'))]
flow <- flow[order(date)][,`:=`(dv_max30d = frollapply(dv,30,max,na.rm=T,align='right'))]


## Load and collate raster stacks of daily MODIS NIR
vec_dates <- fread("data/gee_SAWater_turbidity/DATES_MOD09GA_sur_refl_b02_daily_2000-03-01_2012-12-31.csv")
vec_dates2 <- fread("data/gee_SAWater_turbidity/DATES_MOD09GA_sur_refl_b02_daily_2012-12-31_2023-07-01.csv")

ic <- rast("data/gee_SAWater_turbidity/MOD09GA_sur_refl_b02_daily_2000-03-01_2012-12-31.tif")
ic2 <- rast("data/gee_SAWater_turbidity/MOD09GA_sur_refl_b02_daily_2012-12-31_2023-07-01.tif")

time(ic) <- vec_dates$date
time(ic2) <- vec_dates2$date
names(ic) <- vec_dates$date
names(ic2) <- vec_dates2$date

nir <- c(ic,ic2)
rm(ic2)

## Load and collate raster stacks of daily MODIS RED
vec_dates <- fread("data/gee_SAWater_turbidity/DATES_MOD09GA_sur_refl_b01_daily_2000-03-01_2012-12-31.csv")
vec_dates2 <- fread("data/gee_SAWater_turbidity/DATES_MOD09GA_sur_refl_b01_daily_2012-12-31_2023-07-01.csv")

ic <- rast("data/gee_SAWater_turbidity/MOD09GA_sur_refl_b01_daily_2000-03-01_2012-12-31.tif")
ic2 <- rast("data/gee_SAWater_turbidity/MOD09GA_sur_refl_b01_daily_2012-12-31_2023-07-01.tif")

time(ic) <- vec_dates$date
time(ic2) <- vec_dates2$date
names(ic) <- vec_dates$date
names(ic2) <- vec_dates2$date

red <- c(ic,ic2)
rm(ic2)


## examine LT estimates
ic_u <- mean(ic,na.rm=T)
# ic_50 <- median(ic,na.rm=T)
# ic_sd <- terra::stdev(ic,na.rm=T)

# plot(ic_u, col = inferno(10))
# plot(ic_50, col = inferno(10))
# plot(ic_sd, col = inferno(10))

## Define coords
xy <- ic_u %>% 
  as.data.table(xy=T) %>% 
  select(x,y)

## extract depth and AMC zones
bath <- vect("data/Bathymetry/DetailedModel_DEPTH.shp")
bath <- rasterize(bath, ic, field="Val_1")
names(bath) <- "depth"

amc <- vect("data/AMC/AnalysisExtentMask_GDA94z54.shp")
amc <- terra::project(amc, crs(ic))
amc <- rasterize(amc, ic, field="zone")

vv_bath <- terra::extract(bath, xy[,.(x,y)], xy=T) %>% 
  as.data.table()
vv_amc <- terra::extract(amc, xy[,.(x,y)], xy=T) %>% 
  as.data.table()

tmp_red <- terra::extract(ic, xy)
tmp_red <- cbind(xy,tmp)

dat <- cbind(tmp, vv_bath %>% select(depth),
             vv_amc %>% select(zone))


dat <- dat %>% 
  as_tibble() %>% 
  pivot_longer(-c("x","y","ID","depth","zone"), 
               names_to = c("date"),
               values_to = "red") %>% 
  as.data.table()



library(stars)
cdat <- c(st_as_stars(nir) %>% 
  set_names("nir"),
  st_as_stars(red) %>% 
  set_names("red")) %>% 
  as.data.table(xy=T)

cdat <- merge(cdat, 
      merge(vv_bath,vv_amc,by=c("x","y","ID")),
      by=c("x","y"))

cdat <- cdat[is.na(zone)==F]

cdat[,`:=`(ndvi = (nir-red)/(nir+red))]
cdat[,`:=`(r_nr = red/(nir+red))]
cdat[,`:=`(n_nr = nir/(nir+red))]


cdat <- cdat[,`:=`(g = case_when(
  depth < 10 ~ paste0(str_extract(zone, "\\w*"),"_lt10"),
  depth > 10 ~ paste0(str_extract(zone, "\\w*"),"_gt10")
))][is.na(g)==F]


cdat <- cdat %>% rename(date = time)
cdat[,`:=`(date = date(date))]

sm <- cdat[,.(ndvi50 = median(ndvi,na.rm=T), 
              r_nr50 = median(r_nr,na.rm=T),
              n_nr50 = median(n_nr,na.rm=T), 
              r50  = median(red,na.rm=T),
            n50 = median(nir,na.rm=T)), 
           by=.(date,g)]


sm[,`:=`(date = ymd(date))]
sm[,`:=`(year=year(date),
          month = month(date))]

sm[year==2016][month %in% c(8:12)] %>% 
  ggplot(aes(date, n50, color = g)) +
  geom_line()
cdat[sample(1e5)] %>% ggplot(aes(nir,red))+geom_point()+geom_abline(col='red')



cdat[date %in% seq(ymd("2016-10-01"),ymd("2016-11-01"),by="1 day")] %>% 
  ggplot(aes(x,y,fill=red*scale_factor))+
  geom_tile()+
  coord_sf()+
  scale_fill_viridis_c(option = 'H', 
                       limits = c(0,0.2),
                       oob = scales::squish) + 
  labs(fill = "Red Reflectance", 
       title = "MOD09GA - Terra, Adelaide Metro. Coastal Zones, October 2016", 
       x = NULL,
       y = NULL) +
  facet_wrap(~date, 
             nrow = 3) + 
  theme_linedraw()+
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(size = 7),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = 'grey30'), 
        legend.position = 'bottom')
ggsave(filename = "figures/MOD09GA_2016-10_red.png",
       width = 25,
       height = 20,
       units = 'cm',
       dpi = 350)
