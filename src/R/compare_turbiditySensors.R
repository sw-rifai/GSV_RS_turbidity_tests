pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork)


sn <- fread("data/DEW-Meakin_turbidity/dew_sensor_coordinates.csv")


clean_dew_dat <- function(x){
  tmp <- readxl::read_excel(x, 
                            sheet = "Outliers removed")
  names(tmp) <- tolower(names(tmp))
  tmp <- tmp %>% 
    mutate(time = ymd_hms(paste(mdy(date), time..hh.mm.ss.)))
  return(tmp)
}


flist <- list.files("data/DEW-Meakin_turbidity/B1/",full.names = T, pattern = ".xlsx")


aa <- list()
for(i in 1:length(flist)){
  aa[[i]] <- clean_dew_dat(flist[i])
}

x <- flist[i]
