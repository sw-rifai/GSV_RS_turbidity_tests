pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork)

product_name <- "MOD09GA"
out_var_name <- "blue_to_green" 

# Import ============

## Load and collate raster stacks 
flist <- list.files("data/gee_SAWater_turbidity/",
  pattern = product_name) %>% 
  .[str_detect(.,out_var_name)]
flist_r <- flist[str_detect(flist,".tif")]
flist_dates <- flist[str_detect(flist,"DATES")]



vec_dates <- lapply(flist_dates, FUN = function(x){
  fread(file.path("data/gee_SAWater_turbidity/", x)) %>% 
  pull(date) %>% 
  as_date()
})
l_ic <- lapply(flist_r,FUN= function(x){
  rast(file.path("data/gee_SAWater_turbidity/",x))
})

for(i in 1:length(l_ic)){
  time(l_ic[[i]]) <- vec_dates[[i]]
}

r <- rast(l_ic)

amc <- vect("data/AMC/AnalysisExtentMask_GDA94z54.shp")
amc <- terra::project(amc, crs(r))
amc <- amc %>% select(zone)

r <- (crop(r, amc,mask=T))


fp_out <- paste0(product_name,"_",out_var_name,"_AMC_",
         min(time(r)),"_",max(time(r)),".nc")
fp_out 
writeCDF(r,
  filename=file.path("data/processed_RS/",fp_out),
  varname=out_var_name,
  longname=out_var_name,
  compression=6) 
