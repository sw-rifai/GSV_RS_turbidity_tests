pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork) 

myd <- stars::read_ncdf("data/processed_RS/MYD09GA_blue_to_green_AMC_2002-07-04_2023-10-29.nc", 
                         proxy = F)
mod <- stars::read_ncdf("data/processed_RS/MOD09GA_blue_to_green_AMC_2000-03-01_2023-10-29.nc",
                         proxy = F)

myd <- myd %>% filter(time %in% time(mod))
mod <- mod %>% filter(time %in% time(myd))

names(mod) <- "bg_mod"
names(myd) <- "bg_myd"

# mod <- mod %>% 
#   st_as_stars()
# myd <- myd %>% 
#   st_as_stars()


vec_dates <- time(myd)
vec_dates <- vec_dates[time(myd) %in% time(mod)]

vec_dates <- sample(vec_dates, 500)

mod <- mod %>% filter(time %in% vec_dates)
myd <- myd %>% filter(time %in% vec_dates)

tmp <- c(mod,myd) %>% 
  as.data.table(xy=T) %>% 
  drop_na()

names(tmp) <- c("x","y","time","bg_mod","bg_myd")

ss1 <- tmp %>% 
  drop_na() %>% 
  group_by(x,y) %>% 
  summarize(nobs = n(),
            rho = cor(bg_myd,bg_mod)) %>% 
  ungroup() %>% 
  filter(nobs > 10)
  


## PLOT 
amc <- vect("data/AMC/AnalysisExtentMask_GDA94z54.shp")
amc <- terra::project(amc, crs("EPSG:4326"))
amc <- amc %>% select(zone) %>% 
  st_as_sf()

ss1 %>% 
  ggplot(aes(x,y,fill=rho))+
  geom_raster()+
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
    title="Blue/Green corr. from MYD09GA, MOD09GA") + 
  theme_linedraw()+ 
  theme(panel.grid = element_blank(), 
    panel.background = element_rect(fill='grey'),
    axis.text.x = element_blank())

ggsave("figures/blue_to_green_correlation_MYD09GA_MOD09GA.png",
  width= 10,
  height = 20,
  units='cm',
  dpi =350)


# d %>% 
#   drop_na() %>% 
#   select(red_myd,red_mod) %>% 
#   cor()
# 
# d[id %in% sample(unique(d$id),10)] %>% 
#   drop_na() %>% 
#   ggplot(aes(red_myd,red_mod))+
#   geom_point()+
#   geom_abline(col='red') + 
#   coord_equal() + 
#   facet_wrap(~id, nrow = 1)
# 
# 


# =================================================================
myd <- stars::read_ncdf("data/processed_RS/MYD09GA_blue_to_green_AMC_2002-07-04_2023-10-29.nc", 
                        proxy = F)
mod <- stars::read_ncdf("data/processed_RS/MOD09GA_blue_to_green_AMC_2000-03-01_2023-10-29.nc",
                        proxy = F)

myd <- myd %>% filter(time %in% time(mod))
mod <- mod %>% filter(time %in% time(myd))

names(mod) <- "bg_mod"
names(myd) <- "bg_myd"

# mod <- mod %>% 
#   st_as_stars()
# myd <- myd %>% 
#   st_as_stars()


# vec_dates <- time(myd)
# vec_dates <- vec_dates[time(myd) %in% time(mod)]
# 
# vec_dates <- sample(vec_dates, 500)
# 
# mod <- mod %>% filter(time %in% vec_dates)
# myd <- myd %>% filter(time %in% vec_dates)
# 
tmp <- c(mod,myd) %>% 
  as.data.table(xy=T) %>% 
  drop_na()

names(tmp) <- c("x","y","time","bg_mod","bg_myd")

ss2 <- tmp[,.(mod = median(bg_mod,na.rm=T),
       myd = median(bg_myd,na.rm=T)),by=time]
ss2 %>% 
  pivot_longer(-'time') %>% 
  ggplot(aes(time,value,color=name))+
  geom_smooth(formula = y~s(x,bs='ad'),
              method='gam')
  geom_line()

  
ss3 <- tmp[,.(mod = median(bg_mod,na.rm=T),
                myd = median(bg_myd,na.rm=T)),by=.(x,y)]

ss3[,`:=`(mod_myd = mod-myd)]


ss3$mod %>% summary

ss3 %>% 
  select(-mod_myd) %>% 
  pivot_longer(-c("x","y")) %>% 
  ggplot(aes(x,y,fill=value))+
  geom_raster()+
  geom_sf(data=amc,inherit.aes = F,fill='transparent',
          lwd=1,
          color='black') + 
  geom_sf(data=amc,inherit.aes = F,fill='transparent',
          lwd=0.5,
          color='white') + 
  coord_sf(expand = F)+
  scale_fill_viridis_c(direction = -1,
                       limits=c(0.37,1),
                       oob=scales::squish) +
  labs(x=NULL,y=NULL, 
       fill=expression(paste(frac(Blue,Green))), 
       title="Blue/Green from MYD09GA, MOD09GA") + 
  facet_wrap(~name, nrow=1) +
  theme_linedraw()+ 
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill='grey'),
        axis.text.x = element_blank())

ggsave("figures/blue_to_green_median_MYD09GA_MOD09GA.png",
       width= 10,
       height = 20,
       units='cm',
       dpi =350)


ss3 %>% 
  ggplot(aes(x,y,fill=mod-myd))+
  geom_raster()+
  geom_sf(data=amc,inherit.aes = F,fill='transparent',
          lwd=1,
          color='black') + 
  geom_sf(data=amc,inherit.aes = F,fill='transparent',
          lwd=0.5,
          color='white') + 
  coord_sf(expand = F)+
  scico::scale_fill_scico(palette = 'roma',
                          midpoint = 0) + 
  # scale_fill_viridis_c(direction = -1,
  #                      limits=c(0.37,1),
  #                      oob=scales::squish) +
  labs(x=NULL,y=NULL, 
       fill=expression(paste(frac(Blue,Green))), 
       title="Blue/Green Difference",
       subtitle = "(MOD09GA - MYD09GA)") + 
  theme_linedraw()+ 
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill='grey'),
        axis.text.x = element_blank())

ggsave("figures/blue_to_green_difference_MYD09GA_MOD09GA.png",
       width= 10,
       height = 20,
       units='cm',
       dpi =350)



