# var log10_kd490 = image.expression(
#   'a0 + a1*log10(Blue/Green) + a2*(log10(Blue/Green))**2 + a3*(log10(Blue/Green))**3 + a4*(log10(Blue/Green))**4', {
#     'Blue': image.select('sur_refl_b03'), //
#       'Green': image.select('sur_refl_b04'), //
#       'a0': -0.8813, // sensor specific calib param 
#     'a1': -2.0584, //sensor specific calib param  
#     'a2': 2.5878, //sensor specific calib param  
#     'a3': -3.4885,  //sensor specific calib param 
#     'a4': -1.5061, //sensor specific calib param  
#     'pi': Math.PI,
#     'scale_factor': 1//0.0001 //band info
#   }).rename('log10_mod');
# 
# var kd490 = ee.Image(0.0166).add(ee.Image(10).pow(log10_kd490))
# .rename("kd490_mod")

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


curve(fn_kd490,0.4,1)

data.table(bg = seq(0.15,1,length.out=1000)) %>% 
  mutate(MODIS = fn_kd490_mod(bg),
         L8_OLI = fn_kd490_oli(bg)) %>% 
  pivot_longer(-bg) %>% 
  ggplot(aes(bg,value,color=name))+
  geom_line(lwd=1)+
  scale_y_log10() +
  labs(x='Blue/Green',y='kd490',color="SeaDAS parameter set") + 
  cols4all::scale_color_discrete_c4a_cat(palette = "miscs.okabe") + 
  theme_linedraw()+
  theme(legend.position = 'bottom')
ggsave(filename = "figures/SeaDAS_kd490_blue_to_green_response.png",
       width = 15,
       height = 10,
       units='cm')



