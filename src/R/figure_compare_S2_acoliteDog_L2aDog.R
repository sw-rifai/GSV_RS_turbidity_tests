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



fn_dogliotti2015 <- function(p_w){
  # Dogliotti et al 2015, Remote Sensing of Environment
  # params for 645 nm
  # truncates negative vals and reflectance > 15%

  b0 <- 3.97 # linear scaling hacks
  b1 <- 0.55 # 
  
  A_t <- 228.1
  C <- 0.164
  p_w <- p_w
  p_w <- ifelse(p_w < 0, 0, p_w)
  p_w <- ifelse(p_w > 0.15, 0.15, p_w)
  out <- A_t * (p_w) / 
    (1 - ((p_w)/ C))
  out <- b0 + b1*out
  return(out)
}
fn_dogliotti2015 <- Vectorize(fn_dogliotti2015)
curve(fn_dogliotti2015(x), 0,0.1)

data.table(x=seq(0,0.15,length.out=10000)) %>% 
  mutate(y=fn_dogliotti2015(x)) %>% 
  ggplot(aes(x,y))+
  geom_line()+
  scale_x_log10()+
  scale_y_log10()


# Main ===============================================

## S2 L2A hack dogliotti
a1 <- rast("data/processed_RS/S2-composited_red_AMC_2019-01-02_2023-10-28.nc")
a2 <- a1[[time(a1)=="2019-08-10"]]
# hist(a2)
a3 <- app(a2, fn_dogliotti2015)
plot(a3,range=c(0,50))

## Acolite S2 Dogliotti 2015
tmp1 <- rast("data/S2_acolite/S2B_MSI_2019_08_10_00_46_56_T54HTG_L2W.nc")
tmp2 <- tmp1[['TUR_Dogliotti2015']]
tmp3 <- terra::project(tmp2,a3)

dat <- c(tmp3,a3,a2) %>% 
  as.data.table(xy=T) %>% 
  set_names(c("x","y","aco_dog","l2a_dog","red"))

dat %>% 
  filter(y %between% c(-34.9, -34.6)) %>% 
  select(x,y,aco_dog,l2a_dog) %>% 
  drop_na() %>% 
  sample_n(1000) %>% 
  lm(aco_dog ~ l2a_dog, data=.) %>% summary

# ## RGB
# j <- rast("data/S2_acolite/S2B_MSI_2019_08_10_00_46_56_T54HTG_L2R.nc")
# names(j)
# 
# j2 <- j[[c("rhos_665","rhos_559","rhos_492")]]
# 
# j3 <- terra::project(j2, tmp3)



# Make Figures =====================
# p0 <- ggplot()+
#   geom_spatraster_rgb(data=j3,
#                       max_col_value = 0.18) + 
#   labs(title = "S2 Acolite processed, RGB") +
#   coord_sf(
#     expand = F,
#     ylim = c(-34.9, -34.6)
#   )+
#   scale_x_continuous(breaks = c(138.3,138.4,138.5)) + 
#   theme_map()

p1 <- ggplot()+
  geom_spatraster(data=tmp3,
                  maxcell=1e6)+
  scale_fill_viridis_c(option='H',
                       limits = c(0, 100),
                       oob=scales::squish)+
  labs(fill = "FNU", 
       title="S2 Acolite Dogliotti") +
  coord_sf(
    expand = F,
    ylim = c(-34.9, -34.6)
    )+
  scale_x_continuous(breaks = c(138.3,138.4,138.5)) + 
  theme_map();
p2 <- ggplot()+
  geom_spatraster(data=a3, 
                  maxcell = 1e6)+
  labs(fill="FNU",
       title='S2 L2A Dogliotti') + 
  scale_fill_viridis_c(option='H',
                       limits = c(0, 100),
                       oob=scales::squish)+
  coord_sf(
    expand = F,
    ylim = c(-34.9, -34.6)
  )+
  scale_x_continuous(breaks = c(138.3,138.4,138.5)) + 
  theme_map()


p3 <- dat %>% 
  filter(y %between% c(-34.9, -34.6)) %>% 
  select(x,y,aco_dog,l2a_dog) %>% #pull(l2a_kd) %>% summary
  drop_na() %>% 
  sample_n(1000) %>% 
  ggplot(aes(l2a_dog, aco_dog))+
  geom_point(alpha=0.5,
             shape = 16)+
  geom_abline(col='#cf0000',
              lty=2) + 
  geom_smooth(method='lm') +
  labs(x = "S2 L2A Dogliotti", 
       y = "S2 Acolite Dogliotti") + 
  coord_equal() + 
  theme_map();p3

p_out <- (p3|p1|p2) + 
  plot_layout(guides='collect')+
  plot_annotation(tag_levels = 'a',
                  tag_prefix = '(',
                  tag_suffix = ')')

ggsave(p_out, 
       filename = "figures/figure_compare_S2_acoliteDog_L2aDog.png",
       width = 30,
       height = 12,
       units='cm',
       dpi=350,
       device = grDevices::png)
