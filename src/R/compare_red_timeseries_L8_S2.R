pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork)

var_name <- "red"
  
amc <- vect("data/AMC/AnalysisExtentMask_GDA94z54.shp")
amc <- terra::project(amc, crs("EPSG:4326"))
amc <- amc %>% select(zone) %>% 
  st_as_sf()

bath <- read_sf("data/Bathymetry/DetailedModel_DEPTH.shp")
st_crs(bath) <- st_crs("EPSG:4326")
bath <- bath %>% 
  rename(depth = Val_1) %>% 
  select(depth)

bath <- bath %>% 
  mutate(d_class = if_else(depth < 10, "lt10","gt10")) %>% 
  select(d_class)

lt10 <- bath %>% filter(d_class == "lt10") %>% 
  st_union()
gt10 <- bath %>% filter(d_class == "gt10") %>% 
  st_union()

amc_lt10 <- sf::st_intersection(amc, lt10)
amc_gt10 <- sf::st_intersection(amc,gt10)

amc2 <- bind_rows(
  amc_gt10 %>% 
  mutate(d_class = "gt10"), 
  amc_lt10 %>% 
    mutate(d_class = "lt10"))


amc2 <- amc2 %>% mutate(dz = paste0(str_remove(zone," "),"_",d_class)) %>% 
  select(dz)
amc3 <- vect(amc2)


# Extract from Landsat 8 
l8_r <- terra::rast("data/processed_RS/LC08_red_AMC_2013-08-10_2023-10-09.nc")
system.time(d_l8 <- terra::extract(l8_r, amc3,fun='mean',na.rm=T,bind=T))
d_l8 <- d_l8 %>% as.data.table()
names(d_l8) <- c("dz",time(l8_r) %>% as.character())
d_l8 %>% setDT()
d_l8 <- melt(d_l8, id.vars = 1, variable.name='time', value.name = "red")

ss1 <- d_l8
ss1[,`:=`(time=ymd(time))]

ss1 %>% 
  ggplot(aes(time,red,color=dz))+
  geom_line()


# Extract from S2
s2_r <- terra::rast("data/processed_RS/S2-composited_red_AMC_2019-01-02_2023-10-28.nc")
system.time(d_s2 <- terra::extract(s2_r, amc3, fun='mean',na.rm=T,bind=T))
d_s2 <- d_s2 %>% 
  as.data.table()
names(d_s2) <- c("dz",as.character(time(s2_r)))

d_s2 <- melt(d_s2, id.vars = 1, variable.name='time', value.name = "red")

ss2 <- d_s2[,`:=`(time=ymd(time))]

ss2 %>% 
  ggplot(aes(time,red,color=dz))+
  geom_line()


ss3 <- 
  bind_rows(ss1 %>% 
  mutate(sensor = "L8"),
  ss2 %>% 
  mutate(sensor = "S2"))

arrow::write_parquet(ss3,
                     sink="outputs/AMC_zonal_mean_red_L8_S2.parquet",
                     compression = 'snappy')

ss3 %>% 
  arrange(desc(sensor)) %>% 
  filter(time >= ymd("2019-01-01")) %>% 
  ggplot(aes(time,red,color=sensor))+
  geom_point()+
  geom_smooth(se = F,col='black',lwd=1,aes(group=sensor),
              method='gam',
              formula=y~s(x,bs='ad'))+
  geom_smooth(se = F,lwd=0.6,method='gam',
              formula=y~s(x,bs='ad'))+
  # coord_cartesian(ylim=c(0.3,3),expand=F) + 
  scale_x_date(date_labels = "%Y-%m",
               date_breaks = "6 months") +
  scale_color_discrete_c4a_cat(palette = 'misc.okabe') + 
  labs(x=NULL,
       y=expression(paste(Red))) +
  facet_wrap(~dz,nrow=3) + 
  theme_linedraw() +
  theme(panel.grid.minor = element_blank(), 
        legend.position = 'bottom')

ggsave(filename = "figures/red_compare_timeseries_L8_S2_zones.png",
       width = 30,
       height = 15,
       units='cm',
       dpi=350)


# cols4all::c4a_gui()
# 
# 
# 
# 
# l8 <- stars::read_ncdf("data/processed_RS/LC08_blue_to_green_AMC_2013-08-10_2023-10-09.nc", 
#                        proxy= F)
# st_crs(l8) <- st_crs(4326)
# 
# 
# s2 <- stars::read_ncdf("data/processed_RS/S2-composited_blue_to_green_AMC_2019-01-02_2023-10-28.nc", 
#                         proxy = F)
# st_crs(s2) <- st_crs(4326)
# 
# 
# # d_l8 <- l8 %>% as.data.table(xy=T)
# # names(d_l8) <- c("x","y","time","blue_to_green")
# # d_l8 <- drop_na(d_l8)
# # coords <- d_l8[,.(x,y)] %>% unique()
# # coords <- st_as_sf(coords,coords = c("x","y"))
# # 
# # st_crs(coords) <- st_crs("EPSG:4326")
# # tmp <- sf::st_union(coords, amc2["dz"])
# # 
# # 
# # 
# # dz_s <- st_rasterize(amc2["dz"],template=l8[,,,1])
# # 
# # coords <- dz_s %>% 
# #   as.data.table()
# # coords %>% drop_na
# # names(dz_s)
# 
# 
# system.time(jnk <- stars::st_extract(l8, amc2, FUN = mean, na.rm=T)) # 189
# names(jnk)
# jnk2 <- jnk %>% 
#   as.data.table()
# 
# jnk2 %>% 
#   ggplot(aes(time, blue_to_green)) + 
#   geom_point()
# 
# system.time(jnk <- aggregate(l8, amc2, FUN = mean, na.rm=T)) # 45s
# 
# jnk %>% 
#   st_as_sf() %>% 
#   st_drop_geometry() %>% 
#   mutate(group = 1:nrow(.)) %>% 
#   pivot_longer(-'group',
#                names_to = 'time',
#                values_to = 'bg') %>% 
#   mutate(time = ymd(time)) %>% 
#   ggplot(aes(time,bg,color=factor(group)))+
#   geom_point()+
#   geom_line()
# 
# jnk2 <- jnk %>% 
#   as.data.table(xy=T)
# 
# 
# 
# l8
# s2
# 
# l8 <- l8 %>% 
#   filter(time %in% time(s2))
# 
# s2 <- s2 %>% 
#   filter(time %in% time(l8))
# 
# 
# l8 <- l8 %>% st_as_stars()
# s2 <- s2 %>% st_as_stars()
# 
# names(l8) <- "bg_l8"
# names(s2) <- "bg_s2"
# 
# 
# tmp1 <- l8 %>% as.data.table(xy=T)
# tmp2 <- s2 %>% as.data.table(xy=T)
# 
# tmp3 <- merge(tmp1, tmp2, by=c("x","y","time"))
# 
# tmp4 <- tmp3 %>% drop_na()
# 
# tmp4 %>% select(bg_l8,bg_s2) %>% cor()
# 
# ss1 <- tmp4[,.(rho = cor(bg_l8,bg_s2),
#                l8_p50 = median(bg_l8,na.rm=T),
#                s2_p50 = median(bg_s2,na.rm=T),
#                n = .N),
#             by=.(x,y)]
# 
# ss1[n>3]$rho %>% hist()
# 
# 
# ## PLOT 
# amc <- vect("data/AMC/AnalysisExtentMask_GDA94z54.shp")
# amc <- terra::project(amc, terra::crs("EPSG:4326"))
# amc <- amc %>% select(zone) %>% 
#   st_as_sf()
# 
# ss1[n > 3] %>% 
#   ggplot(aes(x,y,fill=rho))+
#   geom_raster()+
#   geom_sf(data=amc,inherit.aes = F,fill='transparent',
#           lwd=1,
#           color='black') + 
#   geom_sf(data=amc,inherit.aes = F,fill='transparent',
#           lwd=0.5,
#           color='white') + 
#   coord_sf(expand = F)+
#   scico::scale_fill_scico(palette = 'roma',
#                           midpoint=0) + 
#   labs(x=NULL,y=NULL, 
#        fill=expression(paste(rho)), 
#        title="Blue/Green corr. from L8 & S2") + 
#   theme_linedraw()+ 
#   theme(panel.grid = element_blank(), 
#         panel.background = element_rect(fill='grey'),
#         axis.text.x = element_blank())
# 
# ggsave("figures/blue_to_green_correlation_L8_S2.png",
#        width= 10,
#        height = 20,
#        units='cm',
#        dpi =350)
# 
# 
# 
# tmp4[sample(.N,1000)] %>% 
#   ggplot(aes(bg_l8,bg_s2))+
#   geom_point()+
#   geom_abline(col='red')
# 
# 
# ss1[sample(.N,1000)] %>% 
#   ggplot(aes(s2_p50,rho))+
#   geom_point()
# 
# ss1[sample(.N,1000)] %>% 
#   ggplot(aes(l8_p50,rho))+
#   geom_point()
# 
# 
# 
# ss1[is.na(rho)==F][sample(.N,5000)] %>% 
#   ggplot(aes(l8_p50,
#              s2_p50,
#              color=rho))+
#   geom_point(size=1)+
#   geom_abline()+ 
#   coord_equal() +
#   # scale_color_continuous_c4a_div(palette = "brewer.rd_yl_bu")+
#   scale_color_continuous_c4a_div(palette = "scico.roma")+
#   labs(x = "Landsat 8, median Blue/Green reflectance",
#        y = "Sentinel-2, median Blue/Green reflenctance",
#        subtitle = "2019-2023",
#        color = expression(paste(rho))) + 
#   theme_linedraw()+
#   theme(panel.grid = element_blank(),
#         legend.position = c(0.98,0.02),
#         legend.justification = c(0.99,0.01),
#         legend.background = element_rect(fill='transparent'))
# 
# 
# ggsave("figures/blue_to_green_long-termMedianComparison_L8_S2.png",
#        width= 12,
#        height = 10,
#        units='cm',
#        dpi =350)
# 
# 
# 
# 
# ss1[n > 3] %>% 
#   pivot_longer(-c("x","y","rho","n")) %>% 
#   ggplot(aes(x,y,fill=value))+
#   geom_raster()+
#   geom_sf(data=amc,inherit.aes = F,fill='transparent',
#           lwd=1,
#           color='black') + 
#   geom_sf(data=amc,inherit.aes = F,fill='transparent',
#           lwd=0.5,
#           color='white') + 
#   coord_sf(expand = F)+
#   scale_fill_viridis_c(option = 'D', 
#                           # limits = c(0,0.0075),
#                           oob=scales::squish) + 
#   # scico::scale_fill_scico(palette = 'batlow', 
#   #                         limits = c(0,0.01),
#   #                         oob=scales::squish) + 
#   labs(x=NULL,y=NULL, 
#        fill="",
#        title="Blue/Green") +
#   facet_wrap(~name, nrow = 1) + 
#   theme_linedraw()+ 
#   theme(panel.grid = element_blank(), 
#         panel.background = element_rect(fill='grey'),
#         axis.text.x = element_blank(), 
#         legend.position = 'bottom', 
#         legend.key.width = unit(0.075,'npc'))
# ggsave("figures/blue_to_green_median_reflectance_L8_S2.png",
#        width= 20,
#        height = 20,
#        units='cm',
#        dpi =350)
# 
# 
