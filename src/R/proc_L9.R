pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork)
getDTthreads() 

dir_data <- "data/gee_SAWater_turbidity/"
out_var <- "red"

flist <- list.files("data/gee_SAWater_turbidity/",
  pattern = paste0("LC09_"))
flist <- flist[str_detect(flist,out_var)]

flist_dates <- flist[str_detect(flist,"DATES")]
flist_r <- flist[!str_detect(flist,"DATES")]

r <- rast(file.path(dir_data,flist_r))
tmp <- fread(file.path(dir_data,flist_dates))

time(r) <- as_date(tmp$date)

# 
amc <- vect("data/AMC/AnalysisExtentMask_GDA94z54.shp")
amc <- terra::project(amc, crs(r))
amc <- amc %>% select(zone)

r <- crop(r, amc,mask=T)

time(r)

fp_out <- paste0("LC09_",out_var,"_AMC_",
         min(time(r)),"_",max(time(r)),".nc")
writeCDF(r,
  filename=file.path("data/processed_RS/",fp_out),
  varname=out_var,
  longname=out_var,
  compression=6) 


