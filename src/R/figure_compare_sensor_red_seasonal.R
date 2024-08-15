pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork)

tmp <- list.files("outputs/","AMC_zonal_mean_red",full.names =T) %>% 
  lapply(., function(x) arrow::read_parquet(x)) %>% 
  rbindlist()

# cols4all::c4a_gui()
tmp %>% 
  mutate(year=year(time),
         month=lubridate::month(time,label = T)) %>% 
  ggplot(aes(month, red,color=sensor))+
  geom_boxplot(outlier.colour = NA,
               fill='grey70') + 
  coord_cartesian(ylim = c(0,0.1)) + 
  labs(x=NULL,
       y="Red") +
  scale_color_discrete_c4a_cat(palette = "misc.okabe", 
                               reverse = T) +
  # scico::scale_color_scico_d(end = 0.9,
  #                            direction = -1) +
  facet_wrap(~dz) + 
  theme_linedraw()+
  theme(panel.grid.minor = element_blank(), 
        legend.position = 'bottom')
ggsave(filename = "figures/compare_sensor_red_seasonal.png",
       width = 25,
       height = 15,
       units='cm')
