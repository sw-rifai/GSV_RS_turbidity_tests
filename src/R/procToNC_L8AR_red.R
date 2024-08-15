pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork)
getDTthreads() 

# dir_data <- "data/gee_SAWater_turbidity/"
# out_var <- "red"
# 
# flist <- list.files("data/gee_SAWater_turbidity/",
#   pattern = paste0("LC08_"))
# flist <- flist[str_detect(flist,out_var)]
# 
# flist_dates <- flist[str_detect(flist,"DATES")]
# flist_r <- flist[!str_detect(flist,"DATES")]
# 
# r <- rast(file.path(dir_data,flist_r))
# tmp <- fread(file.path(dir_data,flist_dates))
# 
# time(r) <- as_date(tmp$date)
# 
# # 
# amc <- vect("data/AMC/AnalysisExtentMask_GDA94z54.shp")
# amc <- terra::project(amc, crs(r))
# amc <- amc %>% select(zone)
# 
# r <- crop(r, amc,mask=T)
# 
# time(r)
# 
# 
# options =====================================
product <- "Landsat8AquaticReflectance"
out_var_name <- "red"
dir_data <- "/media/sami/drive_a/landsat_marine/AMC/"
flist<- list.files(dir_data)
dir_scratch <- "/home/sami/scratch/"

# # pre-loads ===================================
amc <- sf::read_sf('outputs/AMC_10m_split.gpkg')
amc <- vect(amc)
crs_r <- "PROJCRS[\"WGS 84 / UTM zone 54N\",\n    BASEGEOGCRS[\"WGS 84\",\n        DATUM[\"World Geodetic System 1984\",\n            ELLIPSOID[\"WGS 84\",6378137,298.257223563,\n                LENGTHUNIT[\"metre\",1]]],\n        PRIMEM[\"Greenwich\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n        ID[\"EPSG\",4326]],\n    CONVERSION[\"UTM zone 54N\",\n        METHOD[\"Transverse Mercator\",\n            ID[\"EPSG\",9807]],\n        PARAMETER[\"Latitude of natural origin\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8801]],\n        PARAMETER[\"Longitude of natural origin\",141,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8802]],\n        PARAMETER[\"Scale factor at natural origin\",0.9996,\n            SCALEUNIT[\"unity\",1],\n            ID[\"EPSG\",8805]],\n        PARAMETER[\"False easting\",500000,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8806]],\n        PARAMETER[\"False northing\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8807]]],\n    CS[Cartesian,2],\n        AXIS[\"(E)\",east,\n            ORDER[1],\n            LENGTHUNIT[\"metre\",1]],\n        AXIS[\"(N)\",north,\n            ORDER[2],\n            LENGTHUNIT[\"metre\",1]],\n    USAGE[\n        SCOPE[\"Engineering survey, topographic mapping.\"],\n        AREA[\"Between 138°E and 144°E, northern hemisphere between equator and 84°N, onshore and offshore. Japan. Russian Federation.\"],\n        BBOX[0,138,84,144]],\n    ID[\"EPSG\",32654]]"
# rast(file.path(dir_data,flist[1])) %>% crs()
amc <- project(amc, crs_r)
# amc <- vect("outputs/AMC_10m_split.gpkg")
# amc_zone <- terra::rasterize(amc, rast(file.path(dir_data,flist[1])),
#                              field = "dz")
amc_s <- st_as_sf(amc)


tmp_fps <- list.files(file.path(dir_data,flist[1]),recursive = T)
r_ref <- rast(file.path(dir_data, flist[1], tmp_fps[str_detect(tmp_fps,"AR_BAND4")]))
ext_ref <- ext(amc)

# functions =================

fn_dogliotti <- function(p_w){
  # Dogliotti et al 2015, Remote Sensing of Environment
  # Notes: See ACOLITE for most recent parameter calibrations
  A_t <- 378.46
  C <- 0.19905
  out <- A_t * (p_w) / 
    (1 - ((p_w)/ C))
  return(out)
}


fn_extract <- function(fp){
  scale_factor <- 0.00001
  tmp_fps <- list.files(file.path(dir_data,fp),recursive = T)
  r_red <- rast(file.path(dir_data, fp, tmp_fps[str_detect(tmp_fps,"AR_BAND4")]))
  # values(r_red)[1:100]
  
  r_red <- terra::crop(r_red,  terra::project(amc, r_red), mask=T)
  # values(r_red)[1:100]
  
  r_red <- mask(r_red, ifel(r_red < 0, NA,1))
  # values(r_red)[1:100]
  
  r_red <- terra::clamp(r_red, lower=0, upper=0.5/scale_factor)*scale_factor
  # values(r_red)[1:100]
  
  if(crs(r_red) != crs(r_ref)){
    # high_res0 <- app(high_res, function(x) ifelse(is.na(x), 0, x))
    
    r_red <- terra::project(r_red,r_ref)
    r_red <- terra::crop(r_red, ext_ref)
  }
  
  # r_turb <- app(r_red*scale_factor, fn_dogliotti)
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
  
  time(r_red) <- tmp_time
  names(r_red) <- out_var_name
  return(r_red)
}

# rast(list(
#   fn_extract(flist[1]),
#   fn_extract(flist[121]))) %>% plot()
# jnk <- fn_extract(flist[1])
# values(jnk)[1:100]


r <- lapply(flist, fn_extract)
# rast(r)

r1 <- lapply(r,FUN = function(x){
  if(ext(x) != ext(r[[1]])){
    crop(x,r[[1]],extend=T)
  }else{
    x
  }
})

r2 <- rast(r1)

# 
# rast(list(r[[1]],
#      r[[178]]))
# 
# a <- r[[1]]
# b <- r[[178]]
# 
# crop(b,a,extend=T)
# rast(list(a,crop(b,a,extend=T)))

# r2 <- r %>% 
#   lapply(., function(x) terra::ext(x)==ext(r[[1]]))
# r2 %>% 
#   unlist() %>% 
#   
#   table()


r2 <- r2[[match(sort(time(r2)), time(r2))]]
time(r2) %>% plot

fp_out <- paste0("LC08AR_",out_var_name,"_AMC_",
         min(time(r2)),"_",max(time(r2)),".nc")
fp_out

# writeCDF(r2,
#   filename=file.path("data/processed_RS/",fp_out),
#   varname=out_var_name,
#   longname=out_var_name,
#   compression=6)



r3 <- stars::st_as_stars(r2)
names(r3)
st_get_dimension_values(r3,3) %>% plot
system.time(
  stars::write_mdim(r3,
                  filename=file.path(dir_scratch,
                                     'test.nc'),
                  layer = out_var_name
))
# system.time(stars::write_mdim(r3, filename="data/processed/test.nc"))

system(paste0('gdal_translate -co "COMPRESS=DEFLATE" -co "ZLEVEL=3" ',
       file.path(dir_scratch,
                       'test.nc')," ", file.path("data/processed_RS",fp_out)))

# system(paste0("cdo -f nc4c -z zip_6 setmissval,nan copy ",file.path(dir_scratch,
#                                               'test.nc')," ", file.path("data/processed_RS",fp_out)))
file.remove(file.path(dir_scratch,
                      'test.nc'))

# jnk <- rast(file.path("data/processed_RS",fp_out))
# plot(jnk[[1]],col=viridis(100))
# time(jnk)
# 
# cdo -f nc4c -z zip_6 setmissval,9.96921e+36 copy 
# cdo -f nc4c -z zip_6 -m 9.96921e+36 copy test.nc test3.nc
# 
# gdal_translate -co "COMPRESS=DEFLATE" -co "ZLEVEL=3" test.nc test4.nc
