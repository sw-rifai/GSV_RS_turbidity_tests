pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork)

fn_dogliotti2015 <- function(p_w){
  # Dogliotti et al 2015, Remote Sensing of Environment
  # params for 645 nm
  # truncates negative vals and reflectance > 15%
  A_t <- 228.1
  C <- 0.164
  p_w <- ifelse(p_w < 0, 0, p_w)
  # p_w <- ifelse(p_w > 0.15, 0.15, p_w)
  out <- A_t * (p_w) / 
    (1 - ((p_w)/ C))
  return(out)
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







tmp_uv <- rast("/media/sami/drive_a/landsat_marine/AMC/LC080970842016100502T1-SC20231013070219/LC08_L1TP_097084_20161005_20200906_02_T1_AR_BAND1.tif")*0.00001
tmp_b <- rast("/media/sami/drive_a/landsat_marine/AMC/LC080970842016100502T1-SC20231013070219/LC08_L1TP_097084_20161005_20200906_02_T1_AR_BAND2.tif")*0.00001
tmp_g <- rast("/media/sami/drive_a/landsat_marine/AMC/LC080970842016100502T1-SC20231013070219/LC08_L1TP_097084_20161005_20200906_02_T1_AR_BAND3.tif")*0.00001
tmp_r <- rast("/media/sami/drive_a/landsat_marine/AMC/LC080970842016100502T1-SC20231013070219/LC08_L1TP_097084_20161005_20200906_02_T1_AR_BAND4.tif")*0.00001
amc <- vect("data/AMC/AnalysisExtentMask_GDA94z54.shp")
amc <- project(amc, tmp_uv)
amc_bbox <- ext(amc) %>% 
  vect %>% 
  buffer(20000) %>% 
  ext() %>% 
  vect()


tmp <- c(tmp_b,tmp_g,tmp_r)
names(tmp) <- c("blue","green","red")
tmp <- crop(tmp, amc_bbox, mask=T)
tmp_r %>% plot

amc_bbox

tmp_dog <- app(tmp[['red']], fn_dogliotti2015)

p1 <- ggplot()+
  geom_spatraster_rgb(data = tmp, 
                      maxcell = 1e6,
                      r=3, g=2, b=1, 
                      max_col_value = 0.075) + 
  coord_sf(xlim = c(225895, 293525 - 20000),
           ylim = c(-3937535 + 23000, -3810685 - 5000), 
           expand = F) +
  labs(title = "L8-AR RGB") + 
  theme_map()

p2 <- ggplot()+
  geom_spatraster(data = tmp_dog) + 
  scale_fill_viridis_c(option='B', 
                       limits = c(0,20), 
                       oob =scales::squish) + 
  coord_sf(xlim = c(225895, 293525 - 20000),
           ylim = c(-3937535 + 23000, -3810685 - 5000), 
           expand = F) + 
  labs(fill = NULL,
       title = "2016-10-05") +
  theme_map() +
  theme(legend.position = 'none')

tmp2_r <- rast("/media/sami/drive_a/landsat_marine/AMC/LC080970842016110602T1-SC20231013062951/LC08_L1TP_097084_20161106_20200905_02_T1_AR_BAND4.tif")*0.00001
tmp2_r <- crop(tmp2_r, amc_bbox, mask=T)
tmp2_dog <- app(tmp2_r, fn_dogliotti2015)

p3 <- ggplot()+
  geom_spatraster(data = tmp2_dog) + 
  scale_fill_viridis_c(option='B', 
                       limits = c(0,20), 
                       oob =scales::squish) + 
  coord_sf(xlim = c(225895, 293525 - 20000),
           ylim = c(-3937535 + 23000, -3810685 - 5000), 
           expand = F) + 
  labs(fill = "FNU", 
       title = "2016-11-06") + 
  theme_map() + 
  make_fullsize()

p_out <- (p1|p2|p3)
ggsave(p_out, 
       filename = 'figures/figure_L8AR_2016storm.png',
       width = 30,
       height = 16, 
       units='cm',
       dpi = 350, 
       device = grDevices::png)


# 
# # tmp2_r <- tmp[['red']]/10000
# # tmp2_r[tmp2_r > 1] <- NA_real_
# # tmp2_r[tmp2_r < 0] <- NA_real_
# 
# ggplot()+
#   geom_spatraster(data=tmp_dog)+
#   scale_fill_viridis_c(
#     limits=c(0,30),
#                        oob = scales::squish
#     )
# 
# curve(fn_dogliotti2015(x), 0,0.05)
# 
# tmp2_r[[tmp2_r > 1]]$red <- NA_real_
# tmp2_r %>% plot()
# 
# 
# 
# jnk <- rast("/media/sami/drive_a/landsat_marine/AMC/LC080970842019082702T1-SC20231013052737/LC08_L1TP_097084_20190827_20200826_02_T1_AR_BAND4.tif")*0.00001
# jnk <- crop(jnk,amc_bbox,mask=T)
# plot(jnk, range = c(0, 0.075))
