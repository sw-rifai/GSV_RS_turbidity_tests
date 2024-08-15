pacman::p_load(stars, sf, terra, data.table, tidyverse, lubridate,
               viridis, scico, cols4all, 
               patchwork)
# library(exactextractr)


list.files("outputs/",".parquet")
list.files("outputs/","_blue_to_green_")


tmp1 <- arrow::read_parquet("outputs/Landsat8AquaticReflectance_AMC_blue_to_green.parquet") %>% 
  mutate(sensor = "L8AR")
tmp2 <- arrow::read_parquet("outputs/S2-Level2A_AMC_blue_to_green.parquet") %>% 
  mutate(sensor = "S2-L2A")
tmp3 <- arrow::read_parquet("outputs/MCD43A4_AMC_blue_to_green.parquet") %>% 
  mutate(sensor = "MCD43A4")
tmp4 <- arrow::read_parquet("outputs/MOD09GA_AMC_blue_to_green.parquet") %>% 
  mutate(sensor = "MOD09GA")


dat <- rbindlist(list(tmp1,tmp2,tmp3,tmp4))

dat[,`:=`(count_max = max(count,na.rm=T)),by=.(sensor,dz)]
dat[,`:=`(frac_cover = count/count_max)]



(p1 <- dat[frac_cover >= 0.25] %>% 
  .[time >= ymd("2019-01-01")] %>% 
  .[,`:=`(year=year(time),
          month=month(time))] %>% 
  .[,`:=`(time_ym = ymd(paste(year,month,1)))] %>% 
  .[,.(bg = mean(median,na.rm=T)),
    by=.(sensor, dz,time_ym)] %>% 
  ggplot(aes(time_ym, bg,color=sensor))+
  # geom_point() + 
  geom_line(inherit.aes = F,
            aes(time_ym,bg,group=sensor),
            color='black',
            lwd=1.1)+
  geom_line()+
  scale_color_discrete_c4a_cat(palette = 'misc.okabe') + 
  labs(x=NULL,
       y="Monthly Mean: B/G zonal median",
       color="Product: ") +
  facet_wrap(~dz) +
  theme_linedraw() + 
  theme(legend.position = 'bottom',
  panel.grid.minor = element_blank()))


(p2 <- dat[frac_cover >= 0.1] %>% 
  .[time >= ymd("2019-01-01")] %>% 
  .[time <= ymd("2021-06-01")] %>% 
  .[dz == "CentralZone_lt10"] %>% 
  # .[,`:=`(year=year(time),
  #         month=month(time))] %>% 
  # .[,`:=`(time_ym = ymd(paste(year,month,1)))] %>% 
  # .[,.(q75 = mean(q75,na.rm=T)),
  #   by=.(sensor, dz,time_ym)] %>% 
  ggplot(aes(time, median,color=sensor))+
  geom_point(color='black',size=1.5)+
  geom_point(alpha =0.75, size=1.25) +
  scale_color_discrete_c4a_cat(palette = 'misc.okabe') + 
  labs(x=NULL,
       y="B/G reflectance zonal median",
       title = "Central Zone, < 10 m",
       color = "Product: ") +
  # facet_wrap(~dz) +
  theme_linedraw() + 
  theme(panel.grid.minor = element_blank(), 
        legend.position = 'bottom'))

(p2/p1)
ggsave(filename = "figures/figure_compare_blue_to_green_L8AR_S2_MCD43.png",
       width = 25,
       height = 20,
       units='cm',
       dpi = 350)

