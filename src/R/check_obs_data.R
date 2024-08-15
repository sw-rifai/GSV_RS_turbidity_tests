pacman::p_load(lubridate, sf, terra, tidyverse, data.table,
               mgcv, viridis)


## Bathymetry ==================================
bath <- vect('data/Bathymetry/DetailedModel_DEPTH.shp')
plot(bath)

r <- rast(extent = terra::ext(bath), 
          crs = terra::crs(bath), 
          resolution = 0.001)
r_bath <- rasterize(bath, r, field = 'Val_1')
plot(r_bath,col = inferno(100,direction = -1))

# bathymetry ok!


## Turbidity obs =============================

ns <- fread("data/turbidity_obs_tmp/North_Profiles.csv")
names(ns) <- tolower(names(ns))

ms <- fread("data/turbidity_obs_tmp/Middle_Profiles.csv")
names(ms) <- tolower(names(ms))

ss <- fread("data/turbidity_obs_tmp/Sssouth_Profiles.csv")
names(ss) <- tolower(names(ss))

obs <- rbindlist(list(ns,ms,ss))



oz_poly <- rnaturalearth::ne_countries(country = "Australia",
                                       returnclass = "sf", 
                                       scale = 'large')
amc <- sf::read_sf("data/AMC/AnalysisExtentMask_GDA94z54.shp")
amc <- sf::st_transform(amc, sf::st_crs(oz_poly))
amc <- amc %>% mutate(
  fzone = factor(zone, ordered = T,
                 levels = c("North Zone","Central Zone","South Zone"),
                 labels = c("North Zone","Central Zone","South Zone")))

vec_lims <- st_bbox(st_buffer(amc,10000))

ggplot()+
  geom_sf(data = oz_poly)+ 
  geom_sf(data = amc, aes(fill = fzone)) +
  geom_point(data = obs,  aes(gps_longitude_, 
                        gps_latitude_, 
                        color=turbidity_ntu)) + 
  coord_sf(xlim = vec_lims[c("xmin","xmax")],
           ylim = vec_lims[c("ymin","ymax")], 
           crs = st_crs(amc), 
           ndiscr = 3,
           lims_method = "box") + 
  scale_x_continuous(breaks = seq(
    from = round(vec_lims["xmin"],digits=1), 
    to = round(vec_lims['xmax'],digits=1), 
    by = 0.3
  )) +
  scale_color_viridis_c(option = "H", end = 0.9, 
                        limits = c(0,500),
                        oob=scales::squish) + 
  scale_fill_viridis_d(option = "C", end = 0.9) + 
  labs(fill = "Zone") + 
  theme(panel.background = element_rect(fill='lightblue'))


ss %>% 
  ggplot(aes(gps_longitude_, 
             gps_latitude_, 
             color=turbidity_ntu)) + 
  geom_point()
