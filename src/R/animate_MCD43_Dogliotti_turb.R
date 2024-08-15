pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork)

setwd(here::here())

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


## Load and collate raster stacks of daily MODIS red
vec_dates <- fread("data/gee_SAWater_turbidity/DATES_MCD43A4_Nadir_Reflectance_Band2_daily_2002-02-04_2012-12-31_c2e7c3453c422a1824742911bb6864ed.csv")
vec_dates2 <- fread("data/gee_SAWater_turbidity/DATES_MCD43A4_Nadir_Reflectance_Band2_daily_2012-12-31_2023-07-11_c2e7c3453c422a1824742911bb6864ed.csv")

ic <- rast("data/gee_SAWater_turbidity/MCD43A4_Nadir_Reflectance_Band1_daily_2002-02-04_2012-12-31_61e3bb679d81da01099535c47701ca79.tif")
ic2 <- rast("data/gee_SAWater_turbidity/MCD43A4_Nadir_Reflectance_Band1_daily_2012-12-31_2023-07-11_61e3bb679d81da01099535c47701ca79.tif")

time(ic) <- vec_dates$date
time(ic2) <- vec_dates2$date
names(ic) <- vec_dates$date
names(ic2) <- vec_dates2$date

ic <- c(ic,ic2)
rm(ic2)


amc <- vect("data/AMC/AnalysisExtentMask_GDA94z54.shp")
amc <- terra::project(amc, crs(ic))
amc <- amc %>% select(zone)


fn <- function(r){
  A_t <- 378.46
  p_w <- r
  C <- 0.19905
  scale_factor <- 0.0001
  out <- (A_t * (p_w * scale_factor) / (1 - ((p_w * scale_factor )/ C)))
  return(out)
}

ic_t <- ic %>% 
  crop(., amc,mask=T) %>% 
  fn() 

names(ic_t)

# plot(log1p(ic_t[[3]]),
#      colNA='grey',
#      maxcell=1e6,
#      smooth=F,
#      col=pal,
#      type='continuous',
#      range=c(0,5), 
#      plg = list(title="log1p(FNU)")
# )

vec_dates <- seq(min(time(ic)),max(time(ic)),by='1 month') %>% as.Date()

tmp <- vec_dates %>% 
  lapply(., FUN = function(x){
    out <- clamp(ic_t[[time(ic_t) %between% c(x,x+month(1))]],0,100) %>% max()
    time(out) <- x
    names(out) <- x
    return(out)
  })
tmp <- log1p(rast(tmp))

pal <- inferno(150)
animation::saveGIF(
  animate(
    tmp,
    pause=0.2,
    colNA='grey',
    maxcell=1e6,
    smooth=F,
    col=pal,
    type='continuous',
    range=c(0,5), 
    plg = list(title="log1p(FNU)")  ), 
  movie.name = "mcd43_dogliotti_turb.gif", 
  ani.height = 900,
  ani.width = 600)
