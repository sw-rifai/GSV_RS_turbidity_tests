# fig_AquaOceanColor_kd490TrendMap.R
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

# Main --------------------------------------------------------------------

## Imports =======================
# kd490
kd490 <- terra::rast("data/processed_RS/Aqua_MODIS_kd490.nc")
names(kd490) <- time(kd490)


d1 <- kd490 %>% 
  as.data.table(xy=T)

d1 <- d1 %>% 
  pivot_longer(-c('x','y')) %>% 
  set_names(c("x","y","time","kd490")) %>% 
  as.data.table()
d1[,`:=`(time = as_date(time))]
d1[,`:=`(dyear = decimal_date(time))]
d1[,`:=`(dyear_c = dyear - min(d1$dyear))]
d1[,`:=`(year=year(time),
  month=month(time))]
d1[,`:=`(nobs = sum(is.na(kd490)==F)),by=.(x,y)]
d1[,`:=`(kd490_u = mean(kd490,na.rm=T),
         kd490_sd = mean(kd490,na.rm=T)),by=.(x,y,month)]
d1[,`:=`(kd490_anom = kd490 - kd490_u), by=.(x,y,month)]

## regression =============
s1 <- d1[nobs > 100][,.(betas = list(coef(lm(kd490_anom ~ dyear_c))),
                        kd490_ma = mean(kd490,na.rm=T)),
  by=.(x,y)]
s2 <- s1[,`:=`(b0 = unlist(betas)[1],
            b1 = unlist(betas)[2]),by=.(x,y)]



## prep polys ==============================
oz_coast <- rnaturalearth::ne_coastline(returnclass = 'sf',scale = 'large')
oz_coast <- st_transform(oz_coast, crs = st_crs(myd))
oz_coast <- st_crop(oz_coast, st_bbox(myd))
plot(oz_coast[1])

oz_poly <- rnaturalearth::ne_states(country = 'australia',
  returnclass = 'sf') %>% 
  filter(name == "South Australia") %>% 
  select(name) %>% 
  st_transform(., crs = st_crs(myd)) %>% 
  st_crop(., st_bbox(myd))



s2$kd490_ma %>% hist(100)
s2$kd490_ma %>% summary


## Plots ===========================================
p1 <- s2 %>% 
  ggplot(aes(x,y,fill=kd490_ma))+
  geom_sf(data = oz_poly,
    inherit.aes = F,
    fill='grey30') + 
  geom_raster()+
  geom_sf(data = oz_poly,
    inherit.aes = F,
    color='black') + 
  labs(fill = expression(paste(kd490["MA"]))) + 
  scale_fill_viridis_c(option='B', 
    limits = quantile(s2$kd490_ma,c(0.0,0.975),na.rm=T),
    oob =scales::squish) + 
  # scale_fill_scico(palette = "roma",midpoint=0,direction = -1) + 
  coord_sf(expand = F) + 
  theme_map()+
  make_fullsize();p1


p2 <- s2 %>% 
  ggplot(aes(x,y,fill=b1))+
  geom_sf(data = oz_poly,
    inherit.aes = F,
    fill='grey30') + 
  geom_raster()+
  geom_sf(data = oz_poly,
    inherit.aes = F,
    color='black') + 
  labs(fill = expression(paste(Delta~kd490~yr**-1))) + 
  scale_fill_scico(palette = "roma",midpoint=0,direction = -1) + 
  coord_sf(expand = F) + 
  theme_map()+
  make_fullsize();p2


duration <- d1$dyear %>% range %>% diff
p3 <- s2 %>% 
  # mutate(rr = 100*(b1*duration)/kd490_ma) %>% 
  mutate(rc = 100*b1/kd490_ma) %>% 
  ggplot(aes(x,y,fill=rc))+
  geom_sf(data = oz_poly,
    inherit.aes = F,
    fill='grey30') + 
  geom_raster()+
  geom_sf(data = oz_poly,
    inherit.aes = F,
    color='black') + 
  labs(fill = expression(paste('%'~yr**-1))) + 
  scale_fill_scico(palette = "roma",midpoint=0,direction = -1) + 
  coord_sf(expand = F) + 
  theme_map()+
  make_fullsize();p3


p_out <- (p1|p2|p3) + plot_annotation(tag_levels = 'a',tag_suffix = ')',tag_prefix = '(')&theme(  
  legend.margin = margin(2, 2, 2, 2, "cm")*0.01,
  plot.margin = margin(2, 2, 2, 2, "cm")*0.01)


ggsave(p_out,
  filename = "figures/fig_AquaOceanColor_kd490TrendMap.png", 
  device = grDevices::png,
  width = 14*3, 
  height = 16, 
  units = 'cm', 
  scale = 1,
  dpi = 350) 
