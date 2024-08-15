pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork)

fn_dogliotti2015 <- function(p_w){
  # Dogliotti et al 2015, Remote Sensing of Environment
  # params for 645 nm
  # truncates negative vals and reflectance > 15%
  A_t <- 228.1
  C <- 0.164
  p_w <- ifelse(p_w < 0, 0, p_w)
  # p_w <- ifelse(p_w > 0.15, 0.15, p_w)
  out <- A_t * (p_w) / 
    (1 - ((p_w)/ C))
  return(out)
}

make_fullsize <- function() structure("", class = "fullsizebar")

ggplot_add.fullsizebar <- function(obj, g, name = "fullsizebar") {
  h <- ggplotGrob(g)$heights
  panel <- which(grid::unitType(h) == "null")
  panel_height <- unit(1, "npc") - sum(h[-panel])
  panel_height <- panel_height*0.9
  
  g + 
    guides(fill = guide_colorbar(barheight = panel_height,
                                 title.position = "top")) +
    theme(
      # legend.title = element_text(angle = -90, hjust = 0.5)
      legend.title = element_text(angle = 0, hjust = 0),
    )
}

## Theme ============
theme_map <- function (base_size = 19, base_family = "", base_line_size = base_size/22, 
                       base_rect_size = base_size/22) {
  half_line <- base_size/2
  theme_bw(base_size = base_size, base_family = base_family, 
           base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace% 
    theme(axis.text = element_text(colour = "black", size = rel(0.8)), 
          axis.ticks = element_line(colour = "black", 
                                    linewidth = rel(0.5)),
          axis.ticks.length = unit(-0.1,'cm'),
          legend.text = element_text(size = rel(0.9)),
          panel.border = element_rect(fill = NA, colour = "black", 
                                      linewidth = rel(1)), panel.grid = element_line(colour = "black"), 
          panel.grid.major = element_blank(), #element_line(linewidth = rel(0.1)), 
          panel.grid.minor = element_blank(), #element_line(linewidth = rel(0.05)), 
          strip.background = element_rect(fill = "transparent",
                                          color = 'transparent'), 
          strip.text = element_text(colour = "black", size = rel(0.8), 
                                    margin = margin(0.8 * half_line, 0.8 * half_line, 
                                                    0.8 * half_line, 0.8 * half_line)), 
          legend.position = 'right',
          # legend.position = c(0.5,0.025),
          # legend.justification = c(0.5,0.025),
          legend.direction = 'vertical',
          legend.background =  element_rect(fill='transparent', 
                                            color='transparent'),
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          complete = TRUE)
}


oz_poly <- rnaturalearth::ne_countries(country = "Australia",
                                       returnclass = "sf", 
                                       scale = 'large')

amc <- vect("data/AMC/AnalysisExtentMask_GDA94z54.shp")
amc <- terra::project(amc, crs("EPSG:4326"))
amc <- amc %>% select(zone) %>% 
  st_as_sf()

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


amc2 <- amc2 %>% mutate(dz = paste0(str_remove(zone," "),"_",d_class)) %>% 
  select(dz)
amc3 <- vect(amc2)
amc <- st_as_sf(amc3)

amc <- amc %>% mutate(dz = factor(dz,
                      ordered = T,
                      levels = c("NorthZone_gt10", "NorthZone_lt10", 
                                 "CentralZone_gt10", "CentralZone_lt10",  
                                 "SouthZone_gt10", "SouthZone_lt10"),
                      labels = c("NorthZone_gt10", "NorthZone_lt10", 
                                 "CentralZone_gt10", "CentralZone_lt10",  
                                 "SouthZone_gt10", "SouthZone_lt10")))


oz_poly2 <- st_crop(oz_poly, amc %>% st_buffer(10000))

oz_coast <- rnaturalearth::ne_coastline(returnclass = 'sf',scale = 'large')
oz_coast <- st_crop(oz_coast, amc %>% st_buffer(40000))

ggplot()+
  # geom_sf(data=oz_poly2, 
  #         fill='black',
  #         alpha = 0.5)+
  geom_sf(data=oz_coast, 
          alpha = 0.5)+
  geom_sf(data=amc, 
          aes(fill = dz))+
  scale_fill_brewer(type='qual',
                    palette = 2) + 
  labs(fill = "Zones") + 
  coord_sf(expand = F) + 
  theme_map()+
  theme(panel.background = element_rect(fill='lightblue'))

ggsave(filename = "figures/figure_amc_zones.png",
       width = 20,
       height = 18,
       units='cm',
       dpi = 350)
