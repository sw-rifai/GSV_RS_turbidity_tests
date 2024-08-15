# Install these packages
# install.packages(c("sf","terra","tidyverse"))
library(sf); 
library(terra)
library(tidyverse); 

# Note 1: I'm using two of R's dominant spatial frameworks (sf/stars, and terra)
# They both have useful functions, so it's easier to convert spatial classes
# when needed. 

# Note 2: I didn't use the individual polygons you sent. I could merge them, 
# but relabelling the < 10 and >10 m polygons was taking me too long. 



# This is the AMC zones, delineated at the 10 m depth mark
# see commented code at the end for how it was constructed
amc <- sf::read_sf("outputs/AMC_10m_split.gpkg") %>% 
  vect()



# file paths of rasters
fps_r <- list.files("data/Benthic_Geotiffs/",
           full.names = T, # full path
           pattern = ".tif$") # ends with ".tif"


# function that we will apply to the list of filepaths
extract_class_areas <- function(x){
  tmp <- terra::rast(x) # load raster
  amc2 <- terra::project(amc,crs(tmp)) # project amc vector to whatever CRS the raster is in
  r_zones <- terra::rasterize(amc2, tmp, field='dz') # rasterize amc2 for zonal stats
  
  # stat summary 1
  ss1 <- terra::expanse(tmp, 
                 byValue = T, 
                 unit='km', 
                 zones = r_zones) # calc areas by unique values in raster
  ss1$raster <- names(tmp) # add col name for whatever raster name 
  ss1 <- ss1 %>% 
    rename(area_km2 = area) # rename area w/units
  out <- ss1 %>% as_tibble() # cast data.frame to tibble # not important
  return(out)
}

# applies the function to each filepath
res <- lapply(fps_r, extract_class_areas)
res2 <- bind_rows(res) # binds the results into a dataframe

# write
write_csv(res2, "outputs/boone_zonal_stats.csv")

# sanity check plot
res2 %>% 
  ggplot(aes(area_km2, zone, fill=factor(value)))+
  geom_col() + 
  facet_wrap(~raster)





## Constructing the AMC zones w/10 m depth delineation ###############
# amc <- vect("data/AMC/AnalysisExtentMask_GDA94z54.shp")
# amc <- terra::project(amc, crs("EPSG:4326"))
# amc <- amc %>% 
#   st_as_sf() %>%  
#   select(zone) 
# 
# bath <- read_sf("data/Bathymetry/DetailedModel_DEPTH.shp")
# st_crs(bath) <- st_crs("EPSG:4326")
# bath <- bath %>% 
#   rename(depth = Val_1) %>% 
#   select(depth)
# 
# bath <- bath %>% 
#   mutate(d_class = if_else(depth < 10, "lt10","gt10")) %>% 
#   select(d_class)
# 
# lt10 <- bath %>% filter(d_class == "lt10") %>% 
#   st_union()
# gt10 <- bath %>% filter(d_class == "gt10") %>% 
#   st_union()
# 
# amc_lt10 <- sf::st_intersection(amc, lt10)
# amc_gt10 <- sf::st_intersection(amc,gt10)
# 
# amc2 <- bind_rows(
#   amc_gt10 %>% 
#     mutate(d_class = "gt10"), 
#   amc_lt10 %>% 
#     mutate(d_class = "lt10"))
# 
# amc2 <- amc2 %>% 
#   mutate(dz = paste0(str_remove(zone," "),"_",d_class)) %>% 
#   select(dz)
# 
# sf::write_sf(amc2,  "outputs/AMC_10m_split.shp")
