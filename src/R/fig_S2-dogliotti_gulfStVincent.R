# fig_S2-dogliotti_gulfStVincent.R
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
          # axis.text.x = element_blank(),
          complete = TRUE)
}

fn_dogliotti2015 <- function(p_w){
  # Dogliotti et al 2015, Remote Sensing of Environment
  # params for 645 nm
  # truncates negative vals and reflectance > 15%
  A_t <- 228.1
  C <- 0.164
  p_w <- ifelse(p_w < 0, 0, p_w)
  p_w <- ifelse(p_w > 0.15, 0.15, p_w)
  out <- A_t * (p_w) / 
    (1 - ((p_w)/ C))
  return(out)
}

fn_proc_dogliotti <- function(i){
 i <- mask(i$B4, mask=m_swir, maskvalues=0)  
 dog <- app(i, fn_dogliotti2015, cores = 12)
 return(dog) 
}





## Imports =======================
m_swir <- rast("data/gee_SAWater_turbidity/S2_swir_mask.tif")
# plot(m_swir)

i0 <- terra::rast("data/gee_SAWater_turbidity/S2_median_RGB_2020-11-01_2021-02-28.tif")
i1 <- terra::rast("data/gee_SAWater_turbidity/S2_median_RGB_2021-11-01_2022-02-28.tif")
i2 <- terra::rast("data/gee_SAWater_turbidity/S2_median_RGB_2022-11-01_2023-02-28.tif")
i3 <- terra::rast("data/gee_SAWater_turbidity/S2_median_RGB_2023-11-01_2024-02-28.tif")


j <- list(i0,i1,i2,i3) %>% lapply(., fn_proc_dogliotti)
j1 <- rast(j)
names(j1) <- paste0(2020:2023,"-",2021:2024)




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



p_out <- ggplot()+
  geom_spatraster(data=j1)+
  geom_sf(data=oz_poly,
    fill='grey10',
    inherit.aes = F) + 
  scale_fill_viridis_c(option='B',
    limits=c(0.0,15),
    oob=scales::squish) +
  labs(fill='NTU')+
  coord_sf(expand = F) +
  scale_x_continuous(breaks = c(137.5,138,138.5,139))+
  facet_wrap(~lyr,nrow = 1) + 
  theme_map()+
  make_fullsize()


ggsave(p_out,
  filename = "figures/fig_S2-dogliotti_gulfStVincent.png",
  width = 30,
  height = 12,
  units='cm',
  dpi=350,
  device=grDevices::png)

 