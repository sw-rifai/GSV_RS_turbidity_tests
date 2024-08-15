pacman::p_load(stars, sf, terra, data.table, tidyverse, lubridate,
               viridis,scico)
library(exactextractr)

product <- "Landsat8AquaticReflectance"
out_var_name <- "blue_to_green"
dir_data <- "/media/sami/drive_a/landsat_marine/AMC/"
flist<- list.files(dir_data)

amc <- sf::read_sf('outputs/AMC_10m_split.gpkg')
amc <- vect(amc)
crs_r <- "PROJCRS[\"WGS 84 / UTM zone 54N\",\n    BASEGEOGCRS[\"WGS 84\",\n        DATUM[\"World Geodetic System 1984\",\n            ELLIPSOID[\"WGS 84\",6378137,298.257223563,\n                LENGTHUNIT[\"metre\",1]]],\n        PRIMEM[\"Greenwich\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n        ID[\"EPSG\",4326]],\n    CONVERSION[\"UTM zone 54N\",\n        METHOD[\"Transverse Mercator\",\n            ID[\"EPSG\",9807]],\n        PARAMETER[\"Latitude of natural origin\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8801]],\n        PARAMETER[\"Longitude of natural origin\",141,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8802]],\n        PARAMETER[\"Scale factor at natural origin\",0.9996,\n            SCALEUNIT[\"unity\",1],\n            ID[\"EPSG\",8805]],\n        PARAMETER[\"False easting\",500000,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8806]],\n        PARAMETER[\"False northing\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8807]]],\n    CS[Cartesian,2],\n        AXIS[\"(E)\",east,\n            ORDER[1],\n            LENGTHUNIT[\"metre\",1]],\n        AXIS[\"(N)\",north,\n            ORDER[2],\n            LENGTHUNIT[\"metre\",1]],\n    USAGE[\n        SCOPE[\"Engineering survey, topographic mapping.\"],\n        AREA[\"Between 138°E and 144°E, northern hemisphere between equator and 84°N, onshore and offshore. Japan. Russian Federation.\"],\n        BBOX[0,138,84,144]],\n    ID[\"EPSG\",32654]]"
  # rast(file.path(dir_data,flist[1])) %>% crs()
amc <- project(amc, crs_r)
# amc <- vect("outputs/AMC_10m_split.gpkg")
# amc_zone <- terra::rasterize(amc, rast(file.path(dir_data,flist[1])), 
#                              field = "dz")
amc_s <- st_as_sf(amc)


fn_get_vals <- function(fp){
  tmp_fps <- list.files(file.path(dir_data,fp),recursive = T)
  r_blue <- rast(file.path(dir_data, fp, tmp_fps[str_detect(tmp_fps,"AR_BAND2")]))
  r_green <- rast(file.path(dir_data, fp, tmp_fps[str_detect(tmp_fps,"AR_BAND3")]))
  
  r_bg <- r_blue/r_green
  
  r_bg[r_bg<0] <- NA
  r_bg[r_bg>10] <- NA
  
  tmp_time <- str_extract(tmp_fps[str_detect(tmp_fps,"AR_BAND2")],
                          "_\\d{8}_") %>% 
    str_remove(., "_") %>% 
    str_remove(.,"_") %>% 
    ymd()
  

  tmp2 <- exact_extract(r_bg, amc_s, fun = c('mean', 'median', 'count','stdev','quantile'), 
                        quantiles=c(0.01, 0.05, 0.25,0.75, 0.95,0.99),
                        append_cols = 'dz') %>% setDT()
  scale_factor <- 1
  tmp3 <- cbind(tmp2$dz,tmp2[,2:ncol(tmp2)]*scale_factor)
  names(tmp3) <- names(tmp2)
  tmp3$time <- tmp_time
  tmp3$src_file <- fp
  tmp3$var_name <- out_var_name
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


