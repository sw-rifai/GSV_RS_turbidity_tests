pacman::p_load(stars, sf, terra, data.table, tidyverse, lubridate)
library(exactextractr)

product <- "Landsat8AquaticReflectance"
scale_factor <- 1
out_var_name <- "blue_to_green"
fp <- "data/processed_RS/LC08AR_blue_to_green_AMC_2013-04-06_2023-10-09.nc"
r <- rast(fp)
# dir_data <- "/media/sami/drive_a/landsat_marine/AMC/"
# flist <- list.files(dir_data,
#                     recursive = T, 
#                     pattern = "AR_BAND4")
# flist

amc <- sf::read_sf('outputs/AMC_10m_split.gpkg')
amc <- vect(amc)
crs_r <- r %>% crs()
amc <- vect("outputs/AMC_10m_split.gpkg")
amc <- project(amc, crs_r)
amc_s <- st_as_sf(amc)

vec_time <- time(r) %>% sort()


fn_get_vals <- function(in_time){
  # tmp_time <- str_extract(fp,"_\\d{8}_") %>% 
  #   str_remove(., "_") %>% 
  #   str_remove(.,"_") %>% 
  #   ymd()
  
  tmp1 <- r[[time(r)==in_time]]
  tmp1 <- crop(tmp1,amc,mask=T)
  
  tmp2 <- exact_extract(tmp1, amc_s, fun = c('mean', 'median', 'count','stdev','quantile'), 
                        quantiles=c(0.01, 0.05, 0.25,0.75, 0.95,0.99),
                        append_cols = 'dz') %>% setDT()
  tmp3 <- cbind(tmp2$dz,tmp2[,2:ncol(tmp2)]*scale_factor)
  names(tmp3) <- names(tmp2)
  tmp3$time <- in_time
  # tmp3$src_file <- fp
  return(tmp3)
}


library(future.apply)
plan(multicore, workers=50)
out1 <- future_lapply(vec_time, FUN = fn_get_vals)
out2 <- rbindlist(out1)
plan(sequential)



# dat2 %>% 
#   ggplot(aes(time, median,color=dz))+
#   geom_line()+
#   facet_wrap(~dz,nrow=3)

fp_out <- paste0(product,"_AMC_",out_var_name,".parquet")
fp_out
arrow::write_parquet(out2, 
                     sink=file.path("outputs",fp_out),
                     compression = 'snappy')

