pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork) 


l8 <- stars::read_stars("data/processed_RS/LC08_red_AMC_2013-08-10_2023-10-09.nc", 
                       proxy= T)

ar <- stars::read_stars("data/processed_RS/LC08AR_red_AMC_2013-04-06_2023-10-09.nc", 
                        proxy= T)


l8 <- l8 %>% 
  filter(time %in% time(ar))

ar <- ar %>% 
  filter(time %in% time(l8))


l8 <- l8 %>% st_as_stars()
ar <- ar %>% st_as_stars()

ar2 <- st_warp(ar,l8)

names(l8) <- "red_l8"
names(ar2) <- "red_l8ar"


tmp1 <- l8 %>% as.data.table(xy=T)
tmp2 <- ar2 %>% as.data.table(xy=T)

tmp3 <- merge(tmp1, tmp2, by=c("x","y","time"))

tmp4 <- tmp3 %>% drop_na()

tmp4 %>% select(red_l8,red_l8ar) %>% cor()

ss1 <- tmp4[,.(rho = cor(red_l8,red_l8ar),
               l8_p50 = median(red_l8,na.rm=T),
               l8ar_p50 = median(red_l8ar,na.rm=T),
               n = .N),
            by=.(x,y)]




## PLOT 
amc <- vect("data/AMC/AnalysisExtentMask_GDA94z54.shp")
amc <- terra::project(amc, terra::crs("EPSG:4326"))
amc <- amc %>% select(zone) %>% 
  st_as_sf()

ss1[n > 3] %>% 
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
                          # midpoint=0.5, 
                          limits=c(-1,1)) + 
  labs(x=NULL,y=NULL, 
       fill=expression(paste(rho)), 
       title="Red corr. L8AR & L8LaSRC") + 
  theme_linedraw()+ 
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill='grey'),
        axis.text.x = element_blank())

ggsave("figures/red_correlation_L8AR_L8LaSRC.png",
       width= 10,
       height = 20,
       units='cm',
       dpi =350)



# tmp4[sample(.N,1000)] %>% 
#   ggplot(aes(red_l8,red_s2))+
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


ss1[is.na(rho)==F][sample(.N,10000)] %>% 
  ggplot(aes(l8_p50,
             l8ar_p50,
             color=rho))+
  geom_hline(aes(yintercept=0),col='grey',lty=2)+
  geom_vline(aes(xintercept=0),col='grey',lty=2)+
  geom_point(size=1)+
  geom_abline()+ 
  coord_equal(xlim=c(0,0.1),
                  ylim=c(0,0.1))+
  # scale_color_continuous_c4a_div(palette = "brewer.rd_yl_bu")+
  scale_color_continuous_c4a_div(palette = "scico.roma")+
  labs(x = "L8LaSRC, median red reflectance",
       y = "L8AR, median red reflenctance",
       subtitle = "2013-2023",
       color = expression(paste(rho))) + 
  theme_linedraw()+
  theme(panel.grid = element_blank(),
        legend.position = c(0.1,0.98),
        legend.justification = c(0.01,0.99),
        legend.background = element_rect(fill='transparent'))


ggsave("figures/red_long-termMedianComparison_L8AR_L8LaSRC.png",
       width= 12,
       height = 10,
       units='cm',
       dpi =350)



scico::scico_palette_show()

ss1[n > 3] %>% 
  pivot_longer(-c("x","y","rho","n")) %>% 
  ggplot(aes(x,y,fill=value))+
  geom_tile()+
  geom_sf(data=amc,inherit.aes = F,fill='transparent',
          lwd=1,
          color='black') + 
  geom_sf(data=amc,inherit.aes = F,fill='transparent',
          lwd=0.5,
          color='white') + 
  coord_sf(expand = F)+
  scale_fill_viridis_c(option = 'F', 
                          limits = c(0,0.0075),
                          oob=scales::squish) + 
  # scico::scale_fill_scico(palette = 'batlow', 
  #                         limits = c(0,0.01),
  #                         oob=scales::squish) + 
  labs(x=NULL,y=NULL, 
       fill="reflectance",
       title="Red") +
  facet_wrap(~name, nrow = 1) + 
  theme_linedraw()+ 
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill='grey'),
        axis.text.x = element_blank(), 
        legend.position = 'bottom', 
        legend.key.width = unit(0.1,'npc'))
ggsave("figures/red_median_reflectance_L8AR_L8LaSRC.png",
       width= 20,
       height = 20,
       units='cm',
       dpi =350)

