pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork)
getDTthreads()

# Functions ---------------------------------------------------------------
fn_theme <- function (base_size = 12, base_family = "", base_line_size = base_size/22, 
                      base_rect_size = base_size/22) {
  half_line <- base_size/2
  theme_bw(base_size = base_size, base_family = base_family, 
           base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace% 
    theme(axis.text = element_text(colour = "black", size = rel(0.8)), 
          axis.ticks = element_line(colour = "black", 
                                    linewidth = rel(0.5)),
          axis.ticks.length = unit(-0.1,'cm'),
          panel.border = element_rect(fill = NA, colour = "black", 
                                      linewidth = rel(1)), 
          panel.grid = element_line(colour = "black"), 
          panel.grid.major = element_blank(), #element_line(linewidth = rel(0.1)), 
          panel.grid.minor = element_blank(), #element_line(linewidth = rel(0.05)), 
          strip.background = element_rect(fill = "transparent",
                                          color = 'transparent'), 
          strip.text = element_text(colour = "black", size = rel(0.8), 
                                    margin = margin(0.8 * half_line, 0.8 * half_line, 
                                                    0.8 * half_line, 0.8 * half_line)), 
          # legend.position = c(0.5,0.025),
          # legend.justification = c(0.5,0.025),
          # legend.direction = 'horizontal',
          legend.background =  element_rect(fill='transparent', 
                                            color='transparent'),
          complete = TRUE)
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
    theme(legend.title = element_text(angle = 0, hjust = 0.5,
                                      size = rel(1)))
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

get_mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

## Import and process =====================
# AQUA nasa ocean color: dog
myd <- rast("data/processed_RS/Aqua_MODIS_rrs645.nc")
names(myd) <- time(myd)



r_dog <- terra::app(myd, fn_dogliotti2015, cores = 12)
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

d3[,`:=`(year = year(time),
         month = month(time))]
d4 <- d3[,.(dog = mean(dog,na.rm=T),
            kd = mean(kd,na.rm=T),
            nobs = sum(is.na(kd)==F)),
         by=.(x,y,year,month)]
d4[,`:=`(time = ymd(paste(year,month,1)))]

## most frequent month of max kd490
d5 <- d4 %>% 
  group_by(x,y,year) %>% 
  filter(kd == max(kd,na.rm=T)) %>% 
  ungroup() %>% 
  group_by(x,y) %>% 
  summarize(mmode_kd = get_mode(month)) %>% 
  ungroup()

## most frequent month of max dogliotti
d6 <- d4 %>% 
  group_by(x,y,year) %>% 
  filter(dog == max(dog,na.rm=T)) %>% 
  ungroup() %>% 
  group_by(x,y) %>% 
  summarize(mmode_dog = get_mode(month)) %>% 
  ungroup()

d7 <- merge(d5,d6) %>% as.data.table()


oz_poly <- rnaturalearth::ne_countries(country = "Australia",
                                       returnclass = "sf",
                                       scale = 'large')
oz_poly <- st_crop(oz_poly, ext(myd))



# Plotting ----------------------------------------------------------------


scico::scico_palette_show()
p1 <- d5 %>% ggplot(aes(x,y,fill=mmode_kd))+
  geom_sf(data=oz_poly, inherit.aes = F) + 
  geom_raster()+
  geom_sf(data=oz_poly, inherit.aes = F, 
          fill = 'transparent', 
          color= 'black',
          lwd = 1) + 
  scico::scale_fill_scico(palette = 'romaO',midpoint = 6) +
  coord_sf(expand = F, 
           crs = 4326, 
           label_graticule = "W"
  ) +
  labs(x=NULL,
       y= NULL,
       title = "Month of max Kd490",
       subtitle = "MODIS-AQUA NASA Ocean color (2002 - 2023)",
       fill = 'month') + 
  fn_theme() + 
  make_fullsize()

p2 <- d6 %>% ggplot(aes(x,y,fill=mmode_dog))+
  geom_sf(data=oz_poly, inherit.aes = F) + 
  geom_raster()+
  geom_sf(data=oz_poly, inherit.aes = F, 
          fill = 'transparent', 
          color= 'black',
          lwd = 1) + 
  scico::scale_fill_scico(palette = 'romaO',midpoint = 6) +
  coord_sf(expand = F, 
           crs = 4326, 
           label_graticule = "W"
  ) +
  labs(x=NULL,
       y= NULL,
       title = "Month of max Dogliotti2015",
       subtitle = "MODIS-AQUA NASA Ocean color (2002 - 2023)",
       fill = 'month') + 
  fn_theme() + 
  make_fullsize()

p_out <- (p1|p2) + plot_annotation(tag_levels = 'a',
                                   tag_suffix = ')',
                                   tag_prefix = '(')
ggsave(p_out, 
       filename = "figures/figure_month_max_turbidity_MODISoceanColor.png",
       device = grDevices::png,
       width = 26,
       height = 17, 
       units='cm',
       dpi=350)



