pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork)
getDTthreads()
source("src/R/functions_ggplot_map_theme.R")

# Functions ---------------------------------------------------------------
fn_dogliotti <- function(p_w){
  # Parameters from Dogliotti et al 2015, Remote Sensing of Environment
  A_t <- 228.1
  C <- 0.1641
  out <- A_t * (p_w) / 
    (1 - ((p_w)/ C))
  return(out)
}

## Import and process =====================
# AQUA nasa ocean color: dog
myd <- rast("data/processed_RS/Aqua_MODIS_rrs645.nc")
names(myd) <- time(myd)

oz_poly <- rnaturalearth::ne_countries(country = "Australia",
                                       returnclass = "sf",
                                       scale = 'large')
oz_poly <- st_crop(oz_poly, ext(myd))



r_dog <- terra::app(myd, fn_dogliotti, cores = 12)
d1 <- r_dog %>% 
  as.data.table(xy=T)

d1 <- d1 %>% 
  pivot_longer(-c('x','y')) %>% 
  set_names(c("x","y","time","dog")) %>% 
  as.data.table()
d1[,`:=`(time = as_date(time))]


d3 <- d1[is.na(dog)==F][
  ,`:=`(dyear_start = decimal_date(time)-decimal_date(d1$time %>% min))
]

d5 <- d3[time < min(d3$time) + years(21)][,.(
  coef_dog = list(coef(lm(dog~dyear_start)))
),by=.(x,y)][,`:=`(b0_dog = unlist(coef_dog)[1],
                   b1_dog = unlist(coef_dog)[2]),by=.(x,y)]
d5



# PLOTTING ----------------------------------------------------------------

p2 <- d5 %>% 
  # .[b1_dog %between% quantile(d5$b1_dog,c(0.05,0.95),na.rm=T)] %>% 
  ggplot(aes(x,y,fill=b1_dog))+
  geom_sf(data=oz_poly, inherit.aes = F) + 
  geom_raster()+
  geom_sf(data=oz_poly, inherit.aes = F, 
          fill = 'transparent', 
          color= 'black',
          lwd = 1) + 
  scico::scale_fill_scico(palette = 'roma',
                          midpoint = 0, 
                          direction = -1,
                          # limits = c(-0.0001,0.0001),
                          limits = quantile(d5$b1_dog,c(0.05,0.95),na.rm=T),
                          oob =scales::squish) +
  coord_sf(expand = F, 
           crs = 4326, 
           label_graticule = "W"
  ) +
  labs(x=NULL,
       y= NULL,
       title = "Annual Trend Dogliotti Turbidity",
       subtitle = "MODIS-AQUA NASA Ocean color (2002 - 2023)",
       fill = expression(paste(Delta~NTU~yr**-1))
       ) + 
  fn_theme() + 
  make_fullsize()+
  fn_theme()
p2


p3 <- d5 %>% 
  mutate(rel_change = 100*(21*b1_dog)/b0_dog) %>% 
  # filter(rel_change %between% c(-250,250)) %>% 
  ggplot(aes(x,y,fill= rel_change))+
  geom_sf(data=oz_poly, inherit.aes = F) + 
  geom_raster()+
  geom_sf(data=oz_poly, inherit.aes = F, 
          fill = 'transparent', 
          color= 'black',
          lwd = 1) + 
  scico::scale_fill_scico(palette = 'roma',
                          midpoint = 0, 
                          direction = -1,
                          limits = c(-200, 200),
                          oob =scales::squish) +
  coord_sf(expand = F, 
           crs = 4326, 
           label_graticule = "W"
  ) +
  labs(x=NULL,
       y= NULL,
       title = "Dogliotti Turbidity Relative Change",
       subtitle = "MODIS-AQUA NASA Ocean color (2002 - 2023)",
       fill = "%"
  ) + 
  fn_theme() + 
  make_fullsize()+
  fn_theme()
p3

p_out <- (p2|p3) + plot_annotation(tag_prefix = '(',
                        tag_suffix = ')',
                        tag_levels = 'a')

ggsave(p_out, 
       filename = "figures/plot_AQUA-MODIS-oceancolor_dogliotti_trend.png",
       device = grDevices::png,
       width = 12*2,
       height = 15, 
       units='cm',
       dpi=350)



