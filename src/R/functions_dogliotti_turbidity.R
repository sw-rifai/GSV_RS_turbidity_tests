# var fn_t = function(image) {
#   // Description: Calculate turbidity estimate in Formazin Nephelometric Unit (FNU)
#   // Reference: Dogliotti et al 2015, Remote Sensing of Environment
#   // Input: image with red and SWIR1 bands in reflectane units
#   // Notes: See ACOLITE for most recent parameter calibrations
#   
#   var turbidity = image.expression(
#     '(A_t * (p_w * scale_factor) / (1 - ((p_w * scale_factor )/ C)))', {
#       'p_w': image.select('red'), //red band mid wv_len = 645.5nm
#       'A_t': 378.46, //calib param (ACOLITE)
#       'B_t': 0.33, //calib param (ACOLITE)
#       'C': 0.19905, //calib param (ACOLITE)
#       'pi': Math.PI,
#       'scale_factor': 1//0.0001 //band info
#     }).rename('turb_s2');
#   
#   return(image.addBands(turbidity));
# };

fn_dogliotti <- function(p_w){
  # Nechad 2009? Dogliotti et al 2015, Remote Sensing of Environment
  # Notes: See ACOLITE for most recent parameter calibrations
  # https://www.sciencedirect.com/science/article/pii/S0034425719301014#s0150
  # updated parameter table:
  # https://ars.els-cdn.com/content/image/1-s2.0-S0034425719301014-mmc1.pdf
  A_t <- 378.46
  C <- 0.19905
  out <- A_t * (p_w) / 
        (1 - ((p_w)/ C))
  return(out)
}

fn_dogliotti0 <- function(p_w){
  # Dogliotti et al 2015, Remote Sensing of Environment
  A_t <- 228.1
  B_t <- 0.33
  C <- 0.1641
  out <- (A_t*p_w)/(1 - p_w/C)
  return(out)
}

fn_dogliotti_S2A <- function(p_w){
  # Vanhellemont, Q. (2019). 
  # Adaptation of the dark spectrum fitting atmospheric correction for aquatic 
  # applications of the Landsat and Sentinel-2 archives. 
  # Remote Sensing of Environment, 225, 175-192.
  # updated parameter table:
  # https://ars.els-cdn.com/content/image/1-s2.0-S0034425719301014-mmc1.pdf
  
  # wavelength:  
  A_t <- 268.52
  C <- 0.1725
  
  out <- A_t * (p_w) / 
    (1 - ((p_w)/ C))
  return(out)
}



curve(fn_dogliotti0(x),0,0.1)
curve(fn_dogliotti(x),0,0.1,add=T,col='blue')

data.table(red = seq(0,0.1,length.out=1000)) %>% 
  mutate(v0 = fn_dogliotti0(red),
         v1 = fn_dogliotti(red)) %>% 
  pivot_longer(-red) %>% 
  ggplot(aes(red,value,color=name))+
  geom_line(lwd=1)+
  # scale_y_log10() +
  labs(x='Red reflectance',
       y='Dogliotti Turbidity (FNU)',color="parameter set") + 
  cols4all::scale_color_discrete_c4a_cat(palette = "misc.okabe") + 
  theme_linedraw()+
  theme(legend.position = 'bottom')

ggsave(filename = "figures/dogliotti_turbidity.png",
       width = 15,
       height = 10,
       units='cm')



