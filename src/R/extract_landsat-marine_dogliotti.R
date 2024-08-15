pacman::p_load(stars, sf, terra, data.table, tidyverse, lubridate)
library(exactextractr)

# options =====================================
product <- "Landsat8AquaticReflectance"
out_var_name <- "turb_dog" 
dir_data <- "/media/sami/drive_a/landsat_marine/AMC/"
flist<- list.files(dir_data)


# pre-loads ===================================
amc <- sf::read_sf('outputs/AMC_10m_split.gpkg')
amc <- vect(amc)
crs_r <- "PROJCRS[\"WGS 84 / UTM zone 54N\",\n    BASEGEOGCRS[\"WGS 84\",\n        DATUM[\"World Geodetic System 1984\",\n            ELLIPSOID[\"WGS 84\",6378137,298.257223563,\n                LENGTHUNIT[\"metre\",1]]],\n        PRIMEM[\"Greenwich\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n        ID[\"EPSG\",4326]],\n    CONVERSION[\"UTM zone 54N\",\n        METHOD[\"Transverse Mercator\",\n            ID[\"EPSG\",9807]],\n        PARAMETER[\"Latitude of natural origin\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8801]],\n        PARAMETER[\"Longitude of natural origin\",141,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8802]],\n        PARAMETER[\"Scale factor at natural origin\",0.9996,\n            SCALEUNIT[\"unity\",1],\n            ID[\"EPSG\",8805]],\n        PARAMETER[\"False easting\",500000,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8806]],\n        PARAMETER[\"False northing\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8807]]],\n    CS[Cartesian,2],\n        AXIS[\"(E)\",east,\n            ORDER[1],\n            LENGTHUNIT[\"metre\",1]],\n        AXIS[\"(N)\",north,\n            ORDER[2],\n            LENGTHUNIT[\"metre\",1]],\n    USAGE[\n        SCOPE[\"Engineering survey, topographic mapping.\"],\n        AREA[\"Between 138°E and 144°E, northern hemisphere between equator and 84°N, onshore and offshore. Japan. Russian Federation.\"],\n        BBOX[0,138,84,144]],\n    ID[\"EPSG\",32654]]"
  # rast(file.path(dir_data,flist[1])) %>% crs()
amc <- project(amc, crs_r)
# amc <- vect("outputs/AMC_10m_split.gpkg")
# amc_zone <- terra::rasterize(amc, rast(file.path(dir_data,flist[1])), 
#                              field = "dz")
amc_s <- st_as_sf(amc)





# functions =================
toBinary <- function(i){
  valuesbin = paste0(as.integer(rev(intToBits(i)[1:8])), collapse = "")
  valuesbin = unlist(strsplit(valuesbin, split = ""))
  as.numeric(valuesbin) 
}

fn_dogliotti <- function(p_w){
  # Dogliotti et al 2015, Remote Sensing of Environment
  # Notes: See ACOLITE for most recent parameter calibrations
  A_t <- 378.46
  C <- 0.19905
  out <- A_t * (p_w) / 
    (1 - ((p_w)/ C))
  return(out)
}



r <- rast(ncol=5, nrow=5, vals=c(NA, 1:24))
r2 <- app(r, fn_dogliotti,cores=12)
plot(r2)

a <- app(r, fun = function(i, ff) ff(i), cores=12, ff=toB)
a %>% plot()

gdal_utils(util = "info",
  file.path(dir_data, fp, tmp_fps[str_detect(tmp_fps,"QA_PIXEL")]))


fn_get_vals <- function(fp){
  scale_factor <- 0.00001
  tmp_fps <- list.files(file.path(dir_data,fp),recursive = T)
  r_red <- rast(file.path(dir_data, fp, tmp_fps[str_detect(tmp_fps,"AR_BAND4")]))


  r_red <- terra::crop(r_red,  terra::project(amc, r_red))
  
  r_red <- mask(r_red, ifel(r_red < 0, NA,1))
  r_red <- terra::clamp(r_red, lower=0, upper=0.5/scale_factor)
  # mask(r_red, r_red >= 0, inverse=F) %>% summary
  
  # plot(r_red*scale_factor,col=viridis::viridis(100),range=c(0,0.025))
  
  # AR data is already masked
  # r_qa <- rast(file.path(dir_data, fp, tmp_fps[str_detect(tmp_fps,"QA_PIXEL")])) # NO NEED
  # r_wm <- rast(file.path(dir_data, fp, tmp_fps[str_detect(tmp_fps,"WATER_MASK")]))
  
  r_turb <- app(r_red*scale_factor, fn_dogliotti)
  # plot(r_turb,col=viridis::viridis(100),breaks=seq(0,100,length.out=10))
  # plot(r_turb,col=viridis::inferno(100),range=c(0,10))
  
  # r_qa <- app(r_qa, 
  #     fun = function(x){
  #       toBinary <- function(i){
  #         valuesbin = paste0(as.integer(rev(intToBits(i)[1:8])), collapse = "")
  #         valuesbin = unlist(strsplit(valuesbin, split = ""))
  #         as.numeric(valuesbin)
  #       }
  #       t(sapply(x,toBinary))}, 
  #     cores = 12
  # )
  
  tmp_time <- str_extract(tmp_fps[str_detect(tmp_fps,"AR_BAND2")],
                          "_\\d{8}_") %>% 
    str_remove(., "_") %>% 
    str_remove(.,"_") %>% 
    ymd()
  

  tmp2 <- exact_extract(r_turb, amc_s, fun = c('mean', 'median', 'count','stdev','quantile'), 
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

# dat <- list()
# for(i in 112:length(flist)){
#   print(i)
#   dat[[i]] <- fn_get_vals(flist[i])
# }
# 
# j <- fn_get_vals(flist[112])
# j

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


# which(flist==dat2[median==max(median,na.rm=T)]$src_file)

# library(tidyterra)
# ggplot()+
#   geom_spatraster(data=r_red) +
#   scale_fill_viridis_c(
#     # limits=c(0,50), 
#      oob=scales::squish) + 
#   geom_sf(data=amc_s,inherit.aes = F)
# 
