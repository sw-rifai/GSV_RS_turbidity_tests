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

## S2 L2A hack kd490
a1 <- rast("data/processed_RS/S2-composited_blue_to_green_AMC_2019-01-02_2024-03-11.nc")
a2 <- a1[[time(a1)=="2019-08-10"]]
# hist(a2)
a3 <- app(a2, fn_kd490_oli)
plot(a3)

## Acolite S2 Kd490
tmp1 <- rast("data/S2_acolite/S2B_MSI_2019_08_10_00_46_56_T54HTG_L2W.nc")
tmp2 <- tmp1[['qaa_v6_Kd_490']]
tmp3 <- terra::project(tmp2,a3)

dat <- c(tmp3,a3,a2) %>% 
  as.data.table(xy=T) %>% 
  set_names(c("x","y","qaa_kd","l2a_kd","bg"))

dat %>% 
  filter(y %between% c(-34.9, -34.6)) %>% 
  select(x,y,qaa_kd,l2a_kd) %>% 
  drop_na() %>% 
  sample_n(1000) %>% 
  lm(qaa_kd ~ l2a_kd, data=.) %>% summary

## RGB
j <- rast("data/S2_acolite/S2B_MSI_2019_08_10_00_46_56_T54HTG_L2R.nc")
names(j)

j2 <- j[[c("rhos_665","rhos_559","rhos_492")]]

j3 <- terra::project(j2, tmp3)



# Make Figures =====================
p0 <- ggplot()+
  geom_spatraster_rgb(data=j3,
                      max_col_value = 0.18) + 
  labs(title = "S2 Acolite processed, RGB") +
  coord_sf(
    expand = F,
    ylim = c(-34.9, -34.6)
  )+
  scale_x_continuous(breaks = c(138.3,138.4,138.5)) + 
  theme_map()

p1 <- ggplot()+
  geom_spatraster(data=tmp3)+
  scale_fill_viridis_c(option='H',
                       limits = c(0.5, 5),
                       oob=scales::squish)+
  labs(fill = "Kd490", 
       title="S2 Acolite processed") +
  coord_sf(
    expand = F,
    ylim = c(-34.9, -34.6)
    )+
  scale_x_continuous(breaks = c(138.3,138.4,138.5)) + 
  theme_map();
p2 <- ggplot()+
  geom_spatraster(data=a3, 
                  maxcell = 3e6)+
  labs(fill="Kd490",
       title='S2 L2A processed') + 
  scale_fill_viridis_c(option='H',
                       limits = c(0.5, 5),
                       oob=scales::squish)+
  coord_sf(
    expand = F,
    ylim = c(-34.9, -34.6)
  )+
  scale_x_continuous(breaks = c(138.3,138.4,138.5)) + 
  theme_map()


p3 <- dat %>% 
  filter(y %between% c(-34.9, -34.6)) %>% 
  select(x,y,qaa_kd,l2a_kd) %>% #pull(l2a_kd) %>% summary
  drop_na() %>% 
  sample_n(20000) %>% 
  ggplot(aes(l2a_kd, qaa_kd))+
  geom_point(alpha=0.05,
             shape = 16)+
  geom_abline(col='#cf0000',
              lty=2) + 
  geom_smooth(method='lm') +
  labs(x = "S2 L2A Kd490", 
       y = "S2 Acolite Kd490") + 
  coord_equal() + 
  theme_map()

p_out <- (p0|p3)/(p1|p2) + 
  plot_annotation(tag_levels = 'a',
                  tag_prefix = '(',
                  tag_suffix = ')')

ggsave(p_out, 
       filename = "figures/figure_compare_S2_acoliteKd_L2aKd.png",
       width = 20,
       height = 20,
       units='cm',
       dpi=350,
       device = grDevices::png)
