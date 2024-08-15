pacman::p_load(stars, sf, terra, data.table, tidyverse, lubridate,tidyterra,
               viridis, scico, cols4all, 
               patchwork)

amc <- sf::read_sf("outputs/AMC_10m_split.gpkg") %>% 
  vect()
r <- rast("data/gee_SAWater_turbidity/S2-composited_RGB_2018-12-23.tif")

ext(r)

r <- r %>% 
  crop(., ext( 138.329342322232, 138.517521865378, -34.9, -34.7))

png(filename = "figures/S2-L2A_RGB_2018-12-23.png",
    width = 10,
    height = 13,
    units = 'cm',
    res = 350)
plotRGB(r %>% 
          crop(., ext( 138.329342322232, 138.517521865378, -34.9, -34.7)), 
        stretch = 'lin', 
        maxcell = 5e6,
        smooth = F, colNA='grey', 
        zlim = c(0.01, 0.06), 
        zcol = T)
lines(amc, add=T)
sbar(5, xy="bottomleft", 
     divs=4, 
     cex=.8,
     type='bar',
     lonlat=T,
     ticks=TRUE)
dev.off()

p_rgb <- ggplot()+
  geom_spatraster_rgb(data=r, 
                  maxcell = 1e6, 
                  max_col_value = 0.055) + 
  coord_sf(expand = F) + 
  labs(title = "S2 L2A 2018-12-23") + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_blank())



# Dogliotti turbidity ==============================================
fn_dogliotti <- function(p_w){
  # Dogliotti et al 2015, Remote Sensing of Environment
  # Notes: See ACOLITE for most recent parameter calibrations
  A_t <- 378.46
  C <- 0.19905
  out <- A_t * (p_w) / 
    (1 - ((p_w)/ C))
  return(out)
}

pal_turbo <- turbo(100)
r_turb <- app(r[[1]], fn_dogliotti)


plot(r_turb %>% 
          crop(., ext( 138.329342322232, 138.517521865378, -34.9, -34.7)), 
        maxcell = 0.5e6,
        range = c(0,10),
        col = pal_turbo,
        smooth = F, 
     colNA = 'grey',
     main = "Dogliotti Turbidity (FNU)")


p_turb <- ggplot()+
  geom_spatraster(data=r_turb, 
                  maxcell = 1e6) + 
  scale_fill_viridis_c( option='H',
                        limits = c(0,10),
                        oob=scales::censor) + 
  coord_sf(expand = F) + 
  labs(fill = "FNU", 
       title = "Dogliotti turbidity") + 
  theme(legend.position = 'bottom',
        legend.key.width = unit(1,'cm'), 
        panel.grid = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank())


# kd490 =================================================
fn_kd490_oli <- function(bg){
  # https://oceancolor.gsfc.nasa.gov/resources/atbd/kd/
  vec_params <- c(-0.9054,	-1.5245,	2.2392,	-2.4777,	-1.1099)
  
  a0 <- vec_params[1]
  a1 <- vec_params[2]
  a2 <- vec_params[3]
  a3 <- vec_params[4]
  a4 <- vec_params[5]
  log10_kd490 <- a0 + a1*log10(bg) + a2*(log10(bg))**2 + a3*(log10(bg))**3 + a4*(log10(bg))**4
  kd490 <- 0.0166 + 10**log10_kd490
  return(kd490)
}

fn_kd490_mod <- function(bg){
  # https://oceancolor.gsfc.nasa.gov/resources/atbd/kd/
  a0 <- -0.8813
  a1 <- -2.0584
  a2 <- 2.5878
  a3 <- -3.4885
  a4 <- -1.5061
  log10_kd490 <- a0 + a1*log10(bg) + a2*(log10(bg))**2 + a3*(log10(bg))**3 + a4*(log10(bg))**4
  kd490 <- 0.0166 + 10**log10_kd490
  return(kd490)
}


bg <- (r[[3]]/r[[2]])
r_kd490 <- app(bg, fn_kd490_oli)
r_kd490[r_kd490 > 10] <- NA
hist(r_kd490)

plot(log10(r_kd490),
     maxcell = 0.5e6,
     range = c(0,1),
     col = pal_turbo,
     smooth = F, 
     colNA = 'grey')


p_kd490 <- ggplot()+
  geom_spatraster(data=r_kd490, 
                  maxcell = 1e6) + 
  scale_fill_viridis_c( option='H',
                       trans = 'log10', 
                       limits = c(0.075,2),
                       oob=scales::squish) + 
  coord_sf(expand = F) + 
  labs(fill = expression(paste(m**-1)), 
       title= "Kd490 SeaDAS OLI params.") + 
  theme(legend.position = 'bottom',
        legend.key.width = unit(1,'cm'), 
        panel.grid = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank())


p_out <- p_rgb|p_turb|p_kd490

ggsave(p_out, 
       filename = "figures/figure_S2-L2A_RGB_dogliotti_kd490_vertical.png",
       width = 30,
       height = 15, 
       units='cm', 
       dpi = 350)
