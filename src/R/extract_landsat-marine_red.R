pacman::p_load(stars, sf, terra, data.table, tidyverse, lubridate)
library(exactextractr)

product <- "Landsat8AquaticReflectance"
scale_factor <- 0.00001
out_var_name <- "red"
dir_data <- "/media/sami/drive_a/landsat_marine/AMC/"
flist <- list.files(dir_data,
                    recursive = T, 
                    pattern = "AR_BAND4")
flist

amc <- sf::read_sf('outputs/AMC_10m_split.gpkg')
amc <- vect(amc)
crs_r <- rast(file.path(dir_data,flist[1])) %>% crs()
amc <- project(amc, crs_r)
# amc <- vect("outputs/AMC_10m_split.gpkg")
amc_zone <- terra::rasterize(amc, rast(file.path(dir_data,flist[1])), 
                             field = "dz")
amc_s <- st_as_sf(amc)

plot(amc_zone)

fn_get_vals <- function(fp){
  tmp_time <- str_extract(fp,"_\\d{8}_") %>% 
    str_remove(., "_") %>% 
    str_remove(.,"_") %>% 
    ymd()
  
  tmp1 <- rast(file.path(dir_data,fp))
  tmp1 <- crop(tmp1,amc,mask=T)
  
  tmp2 <- exact_extract(tmp1, amc_s, fun = c('mean', 'median', 'count','stdev','quantile'), 
                        quantiles=c(0.01, 0.05, 0.25,0.75, 0.95,0.99),
                        append_cols = 'dz') %>% setDT()
  tmp3 <- cbind(tmp2$dz,tmp2[,2:ncol(tmp2)]*scale_factor)
  names(tmp3) <- names(tmp2)
  tmp3$time <- tmp_time
  tmp3$src_file <- fp
  rm(tmp_time)
  return(tmp3)
}

dat <- flist %>% 
  lapply(., fn_get_vals)

dat2 <- dat %>% 
  rbindlist()

dat2 %>% 
  ggplot(aes(time, median,color=dz))+
  geom_line()+
  facet_wrap(~dz)

fp_out <- paste0(product,"_AMC_",out_var_name,".parquet")
fp_out
arrow::write_parquet(dat2, 
                     sink=file.path("outputs",fp_out),
                     compression = 'snappy')


