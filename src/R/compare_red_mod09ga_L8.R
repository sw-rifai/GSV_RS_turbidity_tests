pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork)

l8 <- rast("data/processed_RS/LC08_red_AMC_2013-08-10_2023-10-09.nc")
mod <- rast("data/processed_RS/MOD09GA_red_AMC_2000-03-01_2023-07-01.nc")

vec_dates <- time(l8)
vec_dates <- vec_dates[time(l8) %in% time(mod)]

mod <- mod[[time(mod) %in% vec_dates]]
l8 <- l8[[time(l8) %in% vec_dates]]

l8 <- terra::project(l8, mod)
mod <- terra::project(mod,mod) # odd, forces stars to read it later

d_l8 <- st_as_stars(l8) %>% 
  as.data.table(xy=T) %>% 
  set_names(c("x","y","time","red_l8"))

d_mod <- st_as_stars(mod) %>% 
  as.data.table(xy=T) %>% 
  set_names(c("x","y","time","red_mod"))

d <- merge(d_l8,d_mod,by=c("x","y","time"))
d[,`:=`(id = .GRP), by=.(x,y)]

ss1 <- d %>% 
  drop_na() %>% 
  group_by(x,y,id) %>% 
  summarize(nobs = n(),
    rho = cor(red_l8,red_mod)) %>% 
  ungroup() %>% 
  as.data.table()


## PLOT 
amc <- vect("data/AMC/AnalysisExtentMask_GDA94z54.shp")
amc <- terra::project(amc, crs(mod))
amc <- amc %>% select(zone) %>% 
  st_as_sf()

ss1 %>% 
  ggplot(aes(x,y,fill=rho))+
  geom_tile()+
  geom_sf(data=amc,inherit.aes = F,fill='transparent',
    lwd=1,
    color='black') + 
  geom_sf(data=amc,inherit.aes = F,fill='transparent',
    lwd=0.5,
    color='white') + 
  coord_sf(expand = F)+
  scico::scale_fill_scico(palette = 'roma',
    midpoint=0) + 
  labs(x=NULL,y=NULL, 
    fill=expression(paste(rho)), 
    title="Red corr. from L8 and MOD09GA") + 
  theme_linedraw()+ 
  theme(panel.grid = element_blank(), 
    panel.background = element_rect(fill='grey'),
    axis.text.x = element_blank())

ggsave("figures/red_correlation_L8_MOD09GA.png",
  width= 10,
  height = 20,
  units='cm',
  dpi =350)


# d %>% 
#   drop_na() %>% 
#   select(red_l8,red_mod) %>% 
#   cor()
# 
# d[id %in% sample(unique(d$id),10)] %>% 
#   drop_na() %>% 
#   ggplot(aes(red_l8,red_mod))+
#   geom_point()+
#   geom_abline(col='red') + 
#   coord_equal() + 
#   facet_wrap(~id, nrow = 1)
# 
# 
