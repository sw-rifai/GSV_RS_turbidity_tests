pacman::p_load(stars, sf, terra, data.table, tidyverse, lubridate,
               tidyterra,
               mgcv,
               viridis, scico, cols4all, 
               patchwork)


# Functions ---------------------------------------------------------------
theme_map <- function (base_size = 12, base_family = "", base_line_size = base_size/22, 
                       base_rect_size = base_size/22) {
  half_line <- base_size/2
  theme_bw(base_size = base_size, base_family = base_family, 
           base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace% 
    theme(axis.text = element_text(colour = "black", size = rel(0.8)), 
          axis.ticks = element_line(colour = "black", 
                                    linewidth = rel(0.5)),
          axis.ticks.length = unit(-0.1,'cm'),
          panel.border = element_rect(fill = NA, colour = "black", 
                                      linewidth = rel(1)), panel.grid = element_line(colour = "black"), 
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



fn_kd490_oli <- function(bg){
  # hack factors for S2 L2A: 
  hf1 <- -0.3
  hf2 <- 0.073
  hf3 <- 2.622
  
  # https://oceancolor.gsfc.nasa.gov/resources/atbd/kd/
  vec_params <- c(-0.9054,	-1.5245,	2.2392,	-2.4777,	-1.1099)
  bg <- bg+hf1
  a0 <- vec_params[1]
  a1 <- vec_params[2]
  a2 <- vec_params[3]
  a3 <- vec_params[4]
  a4 <- vec_params[5]
  log10_kd490 <- a0 + a1*log10(bg) + a2*(log10(bg))**2 + a3*(log10(bg))**3 + a4*(log10(bg))**4
  kd490 <- hf3*(0.0166 + 10**log10_kd490) + hf2
  return(kd490)
}
fn_kd490_oli <- Vectorize(fn_kd490_oli)
curve(fn_kd490_oli(x), 0.6,1.3)


# Main ===============================================

## S2 L2A predict Kd490 from Acolite Kd490 and B:G
a1 <- rast("data/processed_RS/S2-composited_blue_to_green_AMC_2019-01-02_2024-03-11.nc")
a2 <- a1[[time(a1)=="2019-08-10"]]

## Acolite S2 Kd490
tmp1 <- rast("data/S2_acolite/S2B_MSI_2019_08_10_00_46_56_T54HTG_L2W.nc")
tmp2 <- tmp1[['qaa_v6_Kd_490']]
tmp3 <- terra::project(tmp2,a2)

d1 <- c(tmp3,a2) %>% 
  as.data.table(xy=T) %>% 
  set_names(c("x","y","qaa_kd","l2a_bg"))
d1 <- d1 %>% 
  drop_na()

m_kd <- bam(qaa_kd ~ s(l2a_bg,bs='ts',k=6),
                # method='ML',
                family = Gamma(link='log'),
                select = T,
                discrete = T,
                data = d1 %>% sample_n(100000))
        #         data=d1 %>% 
        # filter(y %between% c(-34.9, -34.6)) %>% 
        # select(x,y,qaa_kd,l2a_bg) %>% 
        # drop_na() %>% 
        # sample_n(50000))
summary(m_kd)
plot(m_kd,rug=T)


# hist(a2)
a3 <- a1[[time(a1)>="2019-08-10"]]
a3 <- a3[[time(a3) <= "2020-08-10"]]
names(a3) <- time(a3)
d2 <- a3 %>% as.data.table(xy=T) %>% 
  pivot_longer(-c('x','y')) %>% 
  set_names(c("x","y","time","l2a_bg")) %>% 
  as.data.table()
head(d2)
d2 <- d2[is.na(l2a_bg)==F]
dim(d2)
d2 <- d2 %>% drop_na()
d2[,`:=`(time = ymd(time))]
system.time(
d3 <- d2 %>% 
  mutate(kd_pred = predict(m_kd,newdata=.,type='response',discrete=T))
)

d4 <- merge(d3, d1[,.(x,y,qaa_kd)], by=c("x","y"))
as.numeric(ymd("2020-01-03") - ymd("2019-08-10"))
d4[,`:=`(tsd = as.numeric(time - ymd("2019-08-10")))]

arrow::write_parquet(d4,
                     sink = "data/processed_RS/S2-L2A_Kd490_dredge_impact_2019-08-10_2020-08-09.parquet",
                     compression = 'snappy')

ss1 <- d4[,.(val = quantile(kd_pred,0.95, na.rm=T)),by=.(time,tsd)]
ss1 %>%
  filter(time <= ymd("2020-03-01")) %>% 
  ggplot(aes(tsd,val))+
  geom_line()

# b0 <- bam(kd_pred ~ te(x,y) + te(x,y,by=tsd), 
#           data=d4[tsd < 60][sample(.N, 1e6)],
#           discrete=T,
#           select=F)
# summary(b0)
# plot(b0, scheme=2)
# 
# b1 <- bam(kd_pred ~ s(tsd,bs='ts') + s(x,y), 
#           data=d4[tsd < 60][sample(.N, 1e6)],
#           discrete=T,
#           select=T)
# summary(b1)
# plot(b1)
# 
# b2 <- bam(I(kd_pred - qaa_kd) ~ 
#             s(qaa_kd, bs='ts',k=5) + 
#             te(x,y,by=tsd, k=c(10, 15)), 
#           data=d4[tsd < 90][sample(.N, 1e6)],
#           discrete=T,
#           select=T)
# summary(b2)
# plot(b2,scheme=2,select=2)
# 
# b3 <- bam(kd_pred ~ s(qaa_kd,by=tsd), 
#           data=d4[tsd < 60][sample(.N, 1e6)],
#           discrete=T,
#           select=T)
# summary(b3)
# plot(b3,scheme=2)
# 
# 
vec_tsd <- unique(d4$tsd)
# d4[time==ymd("2019-8-10")]



oz_coast <- rnaturalearth::ne_coastline(returnclass = 'sf',scale = 'large')
oz_coast <- st_transform(oz_coast, crs = st_crs(i1))
oz_coast <- st_crop(oz_coast, st_bbox(i1))
plot(oz_coast[1])

oz_poly <- rnaturalearth::ne_states(country = 'australia',
                                    returnclass = 'sf') %>% 
  filter(name == "South Australia") %>% 
  select(name) %>% 
  st_transform(., crs = st_crs(a3)) %>% 
  st_crop(., st_bbox(a3))

m_swir <- rast("data/gee_SAWater_turbidity/S2_swir_mask.tif")

# Plotting =============================
p_out <- d4[y>-34.9][tsd<=48] %>% 
  ggplot(aes(x,y,fill=kd_pred-qaa_kd))+
  # geom_spatraster(data = m_swir,
  #                 fill = 'black',
  #                 inherit.aes = F) + 
  geom_raster()+
  # geom_point(data=epicenter,aes(x,y), color='red',
  #            inherit.aes = F) +
  geom_sf(data=oz_poly,
          fill='grey30',
          inherit.aes = F) +
  scale_fill_scico(palette = 'roma',
                   direction = -1,
                   na.value = 'transparent',
                   midpoint=0, 
                   limits  = c(-5,5), 
                   oob = scales::squish) + 
  labs(title = "Days since peak turbidity",
       x = NULL,
       y= NULL,
       fill= expression(paste(Delta~Kd490))) + 
  coord_sf(expand = F, 
           ylim = c(-34.88,-34.59), 
           xlim = c(138.35,138.49)) +
  scale_x_continuous(breaks = c(138.4,138.45))+
  facet_wrap(~tsd,nrow = 2) +
  theme_map() + 
  theme(axis.text.x = element_blank()) + 
  make_fullsize()

ggsave(p_out,
       filename = "figures/figure_dredge_impact_tsd_S2_Kd490.png",
       width = 30,
       height = 30,
       scale = 0.7,
       units='cm',
       dpi=350,
       device=grDevices::png)

# epicenter <- d1[x < 138.45][qaa_kd==max(qaa_kd)][,.(x,y)]
# d1[y>-34.9] %>% ggplot(aes(x,y,fill=qaa_kd))+geom_raster()+
#   geom_point(data=epicenter,aes(x,y), color='red',
#              inherit.aes = F) + 
#   scale_fill_viridis_c()+
#   coord_sf()



# names(a3) <- rep("l2a_bg",dim(a3)[3])
# 
# vec_time <- time(a3)
# fn <- function(x){
#   r <- a3[[time(a3)==x]]
#   names(r) <- "l2a_bg"
#   out <- terra::predict(r, m_kd, na.rm=T, type='response')
#   time(out) <- x
#   return(out)
# }
# a4 <- lapply(vec_time,fn)



# system.time(
# a4 <- app(a2, fun = function(x) mgcv::predict.bam(m_kd,newdata=x,type='response'))
# )
# plot(a3)
# mgcv::predict.bam()



# ## RGB
# j <- rast("data/S2_acolite/S2B_MSI_2019_08_10_00_46_56_T54HTG_L2R.nc")
# names(j)
# 
# j2 <- j[[c("rhos_665","rhos_559","rhos_492")]]
# 
# j3 <- terra::project(j2, tmp3)
# 
