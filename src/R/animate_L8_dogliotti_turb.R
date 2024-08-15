pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork,
               exactextractr)

setwd(here::here())
print(getwd())
# ## extract depth and AMC zones ------------------
ic <- rast("data/gee_SAWater_turbidity/LC08_C02_T1_L2_turb_l8_2013-03-18_2023-09-18.tif")
amc <- vect("data/AMC/AnalysisExtentMask_GDA94z54.shp")
amc <- terra::project(amc, crs(ic))
amc <- amc %>% select(zone)
print("loaded amc")
# 
# bath <- vect("data/Bathymetry/DetailedModel_DEPTH.shp")
# names(bath) <- names(bath) %>% tolower()
# bath <- bath %>% 
#   mutate(depth_class = case_when(val_1 < 5 ~ "< 5",
#                                val_1 >5 & val_1 <10 ~"5 - 10",
#                                val_1 > 10 ~ "10"))
# bath <- aggregate(bath,by='depth_class',dissolve=T)
# bath <- bath %>% select(depth_class)
# bath <- bath %>% filter(is.na(depth_class)==F)
# crs(bath) <- crs(amc)
# 
# bath
# 
# tmp
# aggregate(merge(amc,bath),by=c("zone","depth_class"))
# 
# vect(union(amc,bath))
# vect(intersect(amc,bath))
# vect(intersect(bath,amc))
# 
# union(amc,bath)[[1]]$depth_class
# union(amc,bath)[[2]]$depth_class
# union(amc,bath)[[3]]$depth_class
# union(amc,bath)[[4]]$depth_class
# union(amc,bath)[[5]]$depth_class
# union(amc,bath)[[6]]$depth_class
# union(amc,bath)[[7]]$depth_class

# tmp <- union(amc,bath) %>% 
#   vect %>% 
#   aggregate(., by=c("zone","depth_class"))
# tmp %>% filter(is.na(zone)==F) %>% 
#   filter(is.na(depth_class)==F)
# 
# aggregate(tmp) %>% plot
# 
# plot(tmp)
# names(tmp)
# 
# 
# grid <- rast("data/gee_SAWater_turbidity/LC08_C02_T1_L2_turb_l8_2013-03-18_2023-09-18.tif")[[1]]
# bath <- rasterize(bath, grid, field="depth_class")
# 
# amc <- rasterize(amc,grid,field='zone')
# 
# tmp <- c(bath,amc)
# terra::crosstab(tmp)
# 
# 
# 
# as.data.table(xy=T) %>% 
#   .[is.na(zone)==F]
# tmp %>% 
#   as.polygons()
#   plot()



## Load and collate raster stacks of daily MODIS red
vec_dates <- fread("data/gee_SAWater_turbidity/DATES_LC08_C02_T1_L2_turb_l8_2013-03-18_2023-09-18.csv")
ic <- rast("data/gee_SAWater_turbidity/LC08_C02_T1_L2_turb_l8_2013-03-18_2023-09-18.tif")
time(ic) <- vec_dates$date %>% 
  as.Date()
ic <- crop(ic,amc)
ic <- clamp(ic,lower=0,upper=30)
names(ic)



# plot(log1p(ic[[time(ic)==ymd("2016-10-05")]]) %>% clamp(., lower=0,upper=50),
#      colNA='grey',
#      maxcell=1e6,
#      smooth=T,
#      col=pal,
#      type='continuous',
#      range=c(0,4)
# )

names(ic) <- ymd(str_extract(names(ic),"\\d{8}"))
pal <- inferno(150)
animation::saveGIF(
  animate(log1p(ic),
     colNA='grey',
     maxcell=1e6,
     smooth=T,
     col=pal,
     type='continuous',
     range=c(0,4)
), 
 movie.name = "outputs/L8-LaSRC_dogliotti_turb.gif", 
 ani.height = 900,
 ani.width = 600)

