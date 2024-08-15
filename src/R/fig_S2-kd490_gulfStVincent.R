# fig_S2-kd490_gulfStVincent.R
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

fn_kd490_oli <- function(bg){
  # https://oceancolor.gsfc.nasa.gov/resources/atbd/kd/
  vec_params <- c(-0.9054,	-1.5245,	2.2392,	-2.4777,	-1.1099)
  
  a0 <- vec_params[1]
  a1 <- vec_params[2]
  a2 <- vec_params[3]
  a3 <- vec_params[4]
  a4 <- vec_params[5]
  log10_kd490 <- a0 + a1*log10(bg) + a2*(log10(bg))**2 + a3*(log10(bg))**3 + a4*(log10(bg))**4
  kd490 <- 0.0166 + 10**log10_kd490
  return(kd490)
}

curve(fn_kd490_oli(x),0.3,3)

fn_proc_kd490 <- function(i){
 bg <- mask(i$B2/i$B3, mask=m_swir, maskvalues=0)  
 kd490 <- app(bg, fn_kd490_oli, cores = 12)
 return(kd490) 
}





## Imports =======================
m_swir <- rast("data/gee_SAWater_turbidity/S2_swir_mask.tif")
# plot(m_swir)

i0 <- terra::rast("data/gee_SAWater_turbidity/S2_median_RGB_2020-11-01_2021-02-28.tif")
i1 <- terra::rast("data/gee_SAWater_turbidity/S2_median_RGB_2021-11-01_2022-02-28.tif")
i2 <- terra::rast("data/gee_SAWater_turbidity/S2_median_RGB_2022-11-01_2023-02-28.tif")
i3 <- terra::rast("data/gee_SAWater_turbidity/S2_median_RGB_2023-11-01_2024-02-28.tif")


j <- list(i0,i1,i2,i3) %>% lapply(., fn_proc_kd490)
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
    limits=c(0.025,0.2),
    oob=scales::squish) +
  labs(fill='kd490')+
  coord_sf(expand = F) +
  scale_x_continuous(breaks = c(137.5,138,138.5,139))+
  facet_wrap(~lyr,nrow = 1) + 
  theme_map()+
  make_fullsize()


ggsave(p_out,
  filename = "figures/fig_S2-kd490_gulfStVincent.png",
  width = 30,
  height = 12,
  units='cm',
  dpi=350,
  device=grDevices::png)

