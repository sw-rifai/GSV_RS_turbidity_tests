pacman::p_load(terra,
               tidyverse, 
               data.table, 
               lubridate, 
               sf, 
               tidyterra)

flist <- list.files("/media/sami/drive_a/ocean_color/kd490/Aqua_MODIS/", 
                    pattern = ".nc",
                    full.names = T)

flist


amc <- vect("data/AMC/AnalysisExtentMask_GDA94z54.shp")
amc <- project(amc, crs("EPSG:4326"))

vec_ext <- ext(buffer(amc, width = 50000))

myd <- rast(flist)
myd <- crop(myd, vec_ext, mask=T)

time(myd) <- ymd(str_extract(basename(flist), "\\d{8}"))

writeCDF(myd, 
         filename = "data/processed_RS/Aqua_MODIS_kd490.nc")
