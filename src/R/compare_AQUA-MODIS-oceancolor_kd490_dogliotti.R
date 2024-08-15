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
  # Dogliotti et al 2015, Remote Sensing of Environment
  # Notes: See ACOLITE for most recent parameter calibrations
  A_t <- 378.46
  C <- 0.19905
  out <- A_t * (p_w) / 
    (1 - ((p_w)/ C))
  return(out)
}


## Import and process =====================
# AQUA nasa ocean color: dog
myd <- rast("data/processed_RS/Aqua_MODIS_rrs645.nc")
names(myd) <- time(myd)



r_dog <- terra::app(myd, fn_dogliotti, cores = 12)
d1 <- r_dog %>% 
  as.data.table(xy=T)

d1 <- d1 %>% 
  pivot_longer(-c('x','y')) %>% 
  set_names(c("x","y","time","dog")) %>% 
  as.data.table()
d1[,`:=`(time = as_date(time))]


r_kd <- rast("data/processed_RS/Aqua_MODIS_kd490.nc")
names(r_kd) <- time(r_kd)
d2 <- r_kd %>% 
  as.data.table(xy=T) %>% 
  pivot_longer(-c('x','y')) %>% 
  set_names(c("x","y","time","kd")) %>% 
  as.data.table()
d2[,`:=`(time = as_date(time))]

d3 <- merge(d1, d2, by=c("x","y","time"))


d3 <- d3[is.na(dog)==F][is.na(kd)==F][
  ,`:=`(dyear = decimal_date(time))
]


d4 <- d3[,.(rho = cor(dog,kd)),by=.(x,y)]





oz_poly <- rnaturalearth::ne_countries(country = "Australia",
                                       returnclass = "sf",
                                       scale = 'large')
oz_poly <- st_crop(oz_poly, ext(myd))



# PLOTTING ----------------------------------------------------------------

p1 <- d4 %>% 
  ggplot(aes(x,y,fill=rho))+
  geom_sf(data=oz_poly, inherit.aes = F) + 
  geom_raster()+
  geom_sf(data=oz_poly, inherit.aes = F, 
          fill = 'transparent', 
          color= 'black',
          lwd = 1) + 
  scico::scale_fill_scico(palette = 'roma',midpoint = 0) +
  coord_sf(expand = F, 
           crs = 4326, 
           label_graticule = "W"
  ) +
  labs(x=NULL,
       y= NULL,
       title = "kd490 and Dogliotti correlation",
       subtitle = "MODIS-AQUA NASA Ocean color (2002 - 2023)",
       fill = expression(paste(rho))) + 
  fn_theme() + 
  make_fullsize()+
  fn_theme()

ggsave(p1, 
       filename = "figures/compare_AQUA-MODIS-oceancolor_kd490_dogliotti_rho.png",
       device = grDevices::png,
       width = 12,
       height = 15, 
       units='cm',
       dpi=350)



