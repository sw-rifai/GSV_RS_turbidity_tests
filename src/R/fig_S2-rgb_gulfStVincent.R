# fig_S2-rgb_gulfStVincent.R
pacman::p_load(stars, sf, terra, data.table, tidyverse, lubridate,
               viridis, scico, cols4all, 
               patchwork)

# Functions ---------------------------------------------------------------


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
theme_map <- function (base_size = 15, base_family = "", base_line_size = base_size/22, 
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

# Main --------------------------------------------------------------------

## Imports =======================
i0 <- terra::rast("data/gee_SAWater_turbidity/S2_median_RGB_2020-11-01_2021-02-28.tif")
i1 <- terra::rast("data/gee_SAWater_turbidity/S2_median_RGB_2021-11-01_2022-02-28.tif")
i2 <- terra::rast("data/gee_SAWater_turbidity/S2_median_RGB_2022-11-01_2023-02-28.tif")
i3 <- terra::rast("data/gee_SAWater_turbidity/S2_median_RGB_2023-11-01_2024-02-28.tif")


## prep polys ==============================
oz_coast <- rnaturalearth::ne_coastline(returnclass = 'sf',scale = 'large')
oz_coast <- st_transform(oz_coast, crs = st_crs(i1))
oz_coast <- st_crop(oz_coast, st_bbox(i1))
plot(oz_coast[1])

oz_poly <- rnaturalearth::ne_states(country = 'australia',
  returnclass = 'sf') %>% 
  filter(name == "South Australia") %>% 
  select(name) %>% 
  st_transform(., crs = st_crs(i1)) %>% 
  st_crop(., st_bbox(i1))

# RGB(i1) <- c(3,2,1)
# i1_i <- app(i1*255,floor)
# plotRGB(i1_i, r=3,g=2,b=1,stretch='hist')

p0 <- ggplot()+
  geom_spatraster_rgb(data=i0,
    r=3,g=2,b=1,
    maxcell = 1e6,
    max_col_value = 0.125) + 
  geom_sf(data=oz_poly,
    fill='grey10',
    inherit.aes = F) + 
  labs(title = '2020-11 to 2021-02') + 
  coord_sf(expand=F) + 
 theme_map()  

p1 <- ggplot()+
  geom_spatraster_rgb(data=i1,
    r=3,g=2,b=1,
    maxcell = 1e6,
    max_col_value = 0.125) + 
  geom_sf(data=oz_poly,
    fill='grey10',
    inherit.aes = F) + 
  labs(title = '2021-11 to 2022-02') + 
  coord_sf(expand=F) + 
 theme_map()  


p2 <- ggplot()+
  geom_spatraster_rgb(data=i2,
    r=3,g=2,b=1,
    maxcell = 1e6,
    max_col_value = 0.125) + 
  geom_sf(data=oz_poly,
        fill='grey10',
    inherit.aes = F) + 
  labs(title = '2022-11 to 2023-02') + 
  coord_sf(expand=F) + 
 theme_map()  

p3 <- ggplot()+
  geom_spatraster_rgb(data=i3,
    r=3,g=2,b=1,
    maxcell = 1e6,
    max_col_value = 0.125) + 
  geom_sf(data=oz_poly, 
        fill='grey10',
    inherit.aes = F) + 
  labs(title = '2023-11 to 2024-02') + 
  coord_sf(expand=F) + 
 theme_map()  


p_out <- (p0|p1)/(p2|p3) + plot_annotation(tag_levels = 'a',tag_suffix = ')',tag_prefix = '(')&theme(  
  legend.margin = margin(2, 2, 2, 2, "cm")*0.01,
  plot.margin = margin(2, 2, 2, 2, "cm")*0.01)


ggsave(p_out,
  filename = "figures/fig_S2-rgb_gulfStVincent.png", 
  device = grDevices::png,
  width = 14*2, 
  height = 16*2, 
  units = 'cm', 
  scale = 1,
  dpi = 350) 
