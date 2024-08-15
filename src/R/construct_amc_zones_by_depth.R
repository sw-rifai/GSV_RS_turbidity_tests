

amc <- vect("data/AMC/AnalysisExtentMask_GDA94z54.shp")
amc <- terra::project(amc, crs("EPSG:4326"))
amc <- amc %>% 
  st_as_sf() %>%  
  select(zone) 

bath <- read_sf("data/Bathymetry/DetailedModel_DEPTH.shp")
st_crs(bath) <- st_crs("EPSG:4326")
bath <- bath %>% 
  rename(depth = Val_1) %>% 
  select(depth)

bath <- bath %>% 
  mutate(d_class = if_else(depth < 10, "lt10","gt10")) %>% 
  select(d_class)

lt10 <- bath %>% filter(d_class == "lt10") %>% 
  st_union()
gt10 <- bath %>% filter(d_class == "gt10") %>% 
  st_union()

amc_lt10 <- sf::st_intersection(amc, lt10)
amc_gt10 <- sf::st_intersection(amc,gt10)

amc2 <- bind_rows(
  amc_gt10 %>% 
    mutate(d_class = "gt10"), 
  amc_lt10 %>% 
    mutate(d_class = "lt10"))

amc2 <- amc2 %>% 
  mutate(dz = paste0(str_remove(zone," "),"_",d_class)) %>% 
  select(dz)

sf::write_sf(amc2,  "outputs/AMC_10m_split.shp")

