pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork)
getDTthreads()
source("src/R/functions_ggplot_map_theme.R")

# Functions ---------------------------------------------------------------


## Import and process =====================
# AQUA nasa ocean color: dog
myd <- rast("data/processed_RS/Aqua_MODIS_kd490.nc")
names(myd) <- time(myd)

oz_poly <- rnaturalearth::ne_countries(country = "Australia",
                                       returnclass = "sf",
                                       scale = 'large')
oz_poly <- st_crop(oz_poly, ext(myd))


d1 <- myd %>% 
  as.data.table(xy=T)

d1 <- d1 %>% 
  pivot_longer(-c('x','y')) %>% 
  set_names(c("x","y","time","kd")) %>% 
  as.data.table()
d1[,`:=`(time = as_date(time))]


d3 <- d1[is.na(kd)==F][
  ,`:=`(dyear_start = decimal_date(time)-decimal_date(d1$time %>% min))
]

d5 <- d3[time < min(d3$time) + years(21)][,.(
  nobs = .N,
  coef_kd = list(coef(lm(kd~dyear_start)))
),by=.(x,y)][,`:=`(b0_kd = unlist(coef_kd)[1],
                   b1_kd = unlist(coef_kd)[2]),by=.(x,y)]
d5[]



# PLOTTING ----------------------------------------------------------------

p2 <- d5 %>% 
  .[nobs > 300] %>% 
  # .[b1_kd %between% quantile(d5$b1_kd,c(0.05,0.95),na.rm=T)] %>% 
  ggplot(aes(x,y,fill=b1_kd))+
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
                          limits = quantile(d5$b1_kd,c(0.02,0.98),na.rm=T),
                          oob =scales::squish) +
  coord_sf(expand = F, 
           crs = 4326, 
           label_graticule = "W"
  ) +
  labs(x=NULL,
       y= NULL,
       title = "Annual Trend kd490 Turbidity",
       subtitle = "MODIS-AQUA NASA Ocean color (2002 - 2023)",
       fill = expression(paste(Delta~kd490~yr**-1))
  ) + 
  fn_theme() + 
  make_fullsize()+
  fn_theme()
p2


p3 <- d5 %>% 
  .[nobs > 300] %>% 
  mutate(rel_change = 100*(21*b1_kd)/b0_kd) %>% 
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
                          limits = c(-150, 150),
                          oob =scales::squish) +
  coord_sf(expand = F, 
           crs = 4326, 
           label_graticule = "W"
  ) +
  labs(x=NULL,
       y= NULL,
       title = "kd490 Turbidity Relative Change",
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
       filename = "figures/plot_AQUA-MODIS-oceancolor_kd490_trend.png",
       device = grDevices::png,
       width = 12*2,
       height = 15, 
       units='cm',
       dpi=350)



