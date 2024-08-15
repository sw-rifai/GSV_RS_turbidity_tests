pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork)

scale_factor <- 0.0001;
product_name <- "MOD09GA"
out_var_name <- "red"

# Import ============

## Load and collate raster stacks of daily MODIS SMI
vec_dates <- fread("data/gee_SAWater_turbidity/DATES_MOD09GA_sur_refl_b01_daily_2000-03-01_2012-12-31.csv") %>% 
  pull(date) %>% 
  as_date()
vec_dates2 <- fread("data/gee_SAWater_turbidity/DATES_MOD09GA_sur_refl_b01_daily_2012-12-31_2023-07-01.csv") %>% 
  pull(date) %>% 
  as_date()

ic <- rast("data/gee_SAWater_turbidity/MOD09GA_sur_refl_b01_daily_2000-03-01_2012-12-31.tif")
ic2 <- rast("data/gee_SAWater_turbidity/MOD09GA_sur_refl_b01_daily_2012-12-31_2023-07-01.tif")

time(ic) <- vec_dates
time(ic2) <- vec_dates2

amc <- vect("data/AMC/AnalysisExtentMask_GDA94z54.shp")
amc <- terra::project(amc, crs(ic))
amc <- amc %>% select(zone)

ic <- c(ic,ic2)
r <- (crop(ic, amc,mask=T)*scale_factor)


fp_out <- paste0(product_name,"_",out_var,"_AMC_",
         min(time(r)),"_",max(time(r)),".nc")
fp_out 
writeCDF(r,
  filename=file.path("data/processed_RS/",fp_out),
  varname=out_var,
  longname=out_var,
  compression=6) 
