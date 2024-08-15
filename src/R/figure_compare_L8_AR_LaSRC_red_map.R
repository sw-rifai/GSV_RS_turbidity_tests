pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork, 
               scico)

# Available times
# "2013-08-10", "2013-09-11", "2013-09-27", "2013-11-14", "2013-11-30", 
# "2013-12-16", "2014-02-02", "2014-03-06", "2014-03-22", "2014-06-10", 
# "2014-07-12", "2014-08-13", "2014-08-29", "2014-12-03", "2014-12-19", 
# "2015-01-04", "2015-02-05", "2015-02-21", "2015-04-10", "2015-09-17", 
# "2015-10-19", "2015-12-22", "2016-01-07", "2016-02-08", "2016-03-11", 
# "2016-06-15", "2016-07-17", "2016-08-02", "2016-08-18", "2016-10-05", 
# "2016-12-24", "2017-01-25", "2017-02-10", "2017-02-26", "2017-03-14", 
# "2017-04-15", "2017-06-02", "2017-09-22", "2017-10-08", "2017-11-25", 
# "2017-12-11", "2018-01-28", "2018-02-13", "2018-03-01", "2018-06-05", 
# "2018-06-21", "2018-08-08", "2018-08-24", "2018-09-09", "2018-09-25", 
# "2018-10-11", "2019-02-16", "2019-03-04", "2019-03-20", "2019-04-05", 
# "2019-07-10", "2019-08-27", "2019-09-12", "2019-09-28", "2019-10-30", 
# "2020-01-02", "2020-08-29", "2020-11-01", "2020-11-17", "2020-12-03", 
# "2020-12-19", "2021-03-09", "2021-04-26", "2021-06-13", "2021-06-29", 
# "2021-10-19", "2021-11-04", "2021-12-06", "2021-12-22", "2022-03-12", 
# "2022-03-28", "2022-04-13", "2022-05-15", "2022-07-02", "2022-08-03", 
# "2022-08-19", "2022-09-04", "2022-10-06", "2022-11-07", "2022-12-25", 
# "2023-01-10", "2023-01-26", "2023-02-27", "2023-03-15", "2023-04-16", 
# "2023-07-21", "2023-08-22", "2023-09-23", "2023-10-09"


# OPTIONS ===================
fig_fp_out


# preloads 
amc <- sf::read_sf("outputs/AMC_10m_split.gpkg") %>% vect() 
# %>% 
  # project(., crs(
  #   "PROJCRS[\"WGS 84 / UTM zone 54N\",\n    BASEGEOGCRS[\"WGS 84\",\n        DATUM[\"World Geodetic System 1984\",\n            ELLIPSOID[\"WGS 84\",6378137,298.257223563,\n                LENGTHUNIT[\"metre\",1]]],\n        PRIMEM[\"Greenwich\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n        ID[\"EPSG\",4326]],\n    CONVERSION[\"UTM zone 54N\",\n        METHOD[\"Transverse Mercator\",\n            ID[\"EPSG\",9807]],\n        PARAMETER[\"Latitude of natural origin\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8801]],\n        PARAMETER[\"Longitude of natural origin\",141,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8802]],\n        PARAMETER[\"Scale factor at natural origin\",0.9996,\n            SCALEUNIT[\"unity\",1],\n            ID[\"EPSG\",8805]],\n        PARAMETER[\"False easting\",500000,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8806]],\n        PARAMETER[\"False northing\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8807]]],\n    CS[Cartesian,2],\n        AXIS[\"easting\",east,\n            ORDER[1],\n            LENGTHUNIT[\"metre\",1]],\n        AXIS[\"northing\",north,\n            ORDER[2],\n            LENGTHUNIT[\"metre\",1]],\n    ID[\"EPSG\",32654]]"
  # ))


lasrc <- rast("data/processed_RS/LC08_red_AMC_2013-08-10_2023-10-09.nc")
ar <- rast("data/processed_RS/LC08AR_red_AMC_2013-04-06_2023-10-09.nc")

vec_time1 <- time(lasrc)
vec_time2 <- time(ar)


vec_time <- vec_time1[vec_time1 %in% vec_time2]



# pal <- viridis::rocket(100)
pal <- scico(palette = 'lajolla',n = 500,direction = -1)

dir.create("figures/figure_compare_L8_AR_LaSRC_red")

fn_make_fig <- function(target_time){
  lasrc <- rast("data/processed_RS/LC08_red_AMC_2013-08-10_2023-10-09.nc")
  lasrc <- lasrc[[time(lasrc)==target_time]]
  
  
  ar <- rast("data/processed_RS/LC08AR_red_AMC_2013-04-06_2023-10-09.nc")
  ar <- ar[[time(ar)==target_time]]
  ar <- ar %>% project(., crs(lasrc))
  
  ar[ar==0] <- NA
  
  
  fig_fp_out <- paste0("figures/figure_compare_L8_AR_LaSRC_red/figure_compare_L8_AR_LaSRC_red_",target_time,".png")
  
  lim_upper <- global(ar,fun=function(x) quantile(x, 0.995,na.rm=T)) %>% as.numeric()
  # dev.new(width = 10, height = 10, noRStudioGD = TRUE)
  png(filename = fig_fp_out, 
      width = 15, 
      height = 15, 
      units = 'cm',
      res=350)
  par(mfrow=c(1,2))
  plot(ar,range=c(0,lim_upper),col=pal,main="L8 Aquatic Ref.: Red", 
       background = "grey40", 
       mar=c(1.5,1,1.5,4))
  lines(amc)
  plot(lasrc,range=c(0,lim_upper),col=pal,main="L8 LaSRC: Red", 
       background = "grey40", 
       mar=c(1.5,1,1.5,4))
  lines(amc)
  dev.off()
  
}

vec_time %>% 
  lapply(fn_make_fig)
