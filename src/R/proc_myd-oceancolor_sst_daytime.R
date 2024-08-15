pacman::p_load(terra,
               tidyverse, 
               data.table, 
               lubridate, 
               sf, 
               tidyterra) 

flist <- list.files("/media/sami/drive_a/ocean_color/sst_daytime/MODIS_AQUA/", 
                    pattern = ".nc",
                    full.names = T)

flist


amc <- vect("data/AMC/AnalysisExtentMask_GDA94z54.shp")
amc <- project(amc, crs("EPSG:4326"))

vec_ext <- ext(buffer(amc, width = 50000))

myd <- rast(flist)
myd <- crop(myd, vec_ext, mask=T)
myd_sst <- myd[[names(myd)=="sst"]]
myd_qsst <- myd[[names(myd)=="qual_sst"]]

time(myd_sst) <- ymd(str_extract(basename(flist), "\\d{8}"))
time(myd_qsst) <- ymd(str_extract(basename(flist), "\\d{8}"))

writeCDF(myd_sst, 
         filename = "data/processed_RS/Aqua_MODIS_sst_daytime.nc")
writeCDF(myd_qsst, 
         filename = "data/processed_RS/Aqua_MODIS_qual-sst_daytime.nc")
