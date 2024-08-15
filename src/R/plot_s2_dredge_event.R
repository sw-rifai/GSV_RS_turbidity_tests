pacman::p_load(stars, sf, terra, data.table, tidyverse, lubridate,
               arrow, mgcv, gratia, mgcViz, patchwork, 
               tidyterra)

# Functions =========
fn_plot <- function(x){
  plot(x, col=viridis::inferno(20), 
       # breaks=seq(from=0, to=0.05,length.out=20)
       range=c(0,0.05))
}

dt <- as.data.table

# Main =========
d_s2 <- terra::rast("data/processed_RS/S2-composited_red_AMC_2019-01-02_2023-10-28.nc")
d_s2 <- terra::crop(d_s2, terra::ext(c(138.3,138.5,-34.9,-34.75)))

d3 <- d_s2[[time(d_s2) >= ymd("2019-05-15")]]
d3 <- d3[[time(d3) <= ymd("2019-10-29")]]

d_count <- terra::global(d3, fun = "notNA") %>% dt()
d_count[,`:=`(time = time(d3))]

d_keep <- d_count[notNA > max(d_count$notNA)*0.5]

d4 <- d3[[time(d3) %in% d_keep$time]]

names(d4) <- time(d4)
p_out <- ggplot()+
  geom_spatraster(data=d4, 
                  maxcell = 0.5e6) + 
  scale_fill_viridis_c(option='B',
                       limits=c(0,0.1),
                       oob=scales::squish) +
  labs(fill = "Red reflectance") + 
  facet_wrap(~lyr,nrow = 4) + 
  theme_linedraw()+
  theme(legend.position = 'bottom',
        legend.key.width = unit(0.1, 'npc'),
        axis.text.x = element_text(angle=90))

ggsave(p_out,
       filename = "figures/figure_s2_2019dredgePeriod_mapsRed.png",
       width = 30,
       height = 18,
       units='cm',
       dpi=350)

# plot(d4, 
#      col=viridis::inferno(20), 
#      # breaks=seq(from=0, to=0.05,length.out=20)
#      range=c(0,0.05))
# 
# 
# 
# 
# dim(d3)[3]
# 
# 
# 
# tmp1 <- d3[[1:50]] %>% mean(na.rm=T)
# tmp2 <- d3[[51:100]] %>% mean(na.rm=T)
# tmp3 <- d3[[101:150]] %>% mean(na.rm=T)
# 
# fn_plot(tmp1)
# fn_plot(tmp2)
# fn_plot(tmp3)
# 
# 
# d3[[1:10]] %>% mean(na.rm=T) %>% fn_plot() # distinct line
# d3[[11:20]] %>% mean(na.rm=T) %>% fn_plot() # broadly present
# d3[[21:30]] %>% mean(na.rm=T) %>% fn_plot() # present
# d3[[31:40]] %>% mean(na.rm=T) %>% fn_plot() # present, atmos contam?
# d3[[41:50]] %>% mean(na.rm=T) %>% fn_plot() # very bright
# d3[[51:60]] %>% mean(na.rm=T) %>% fn_plot() # bright
# d3[[61:70]] %>% mean(na.rm=T) %>% fn_plot() # attenuating
# 
# time(d3[[1:10]])
# time(d3)
# 

