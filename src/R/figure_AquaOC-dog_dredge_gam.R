pacman::p_load(signal, stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork)
getDTthreads()

# Functions ================================
## Theme ============
fn_theme <- function (base_size = 12, base_family = "", base_line_size = base_size/22, 
                      base_rect_size = base_size/22){
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
          legend.position = c(0.5,0.025),
          legend.justification = c(0.5,0.025),
          legend.direction = 'horizontal',
          legend.background =  element_rect(fill='transparent', 
                                            color='transparent'),
          complete = TRUE)
}

fn_dogliotti2015 <- function(p_w){
  # Dogliotti et al 2015, Remote Sensing of Environment
  # params for 645 nm
  # truncates negative vals and reflectance > 15%
  A_t <- 228.1
  C <- 0.164
  p_w <- ifelse(p_w < 0, 0, p_w)
  p_w <- ifelse(p_w > 0.15, 0.15, p_w)
  out <- A_t * (p_w) / 
    (1 - ((p_w)/ C))
  return(out)
}


# Main ========================================

## Import AQUA ocean color ===========
myd <- terra::rast("data/processed_RS/Aqua_MODIS_rrs645.nc")
names(myd) <- time(myd)
d1 <- as.data.table(myd,xy=T)

d1 <- d1 %>% 
  pivot_longer(-c('x','y')) %>% 
  set_names(c("x","y","time","red")) %>% 
  as.data.table()
d1[,`:=`(time = ymd(time))]
d1[,`:=`(dog = fn_dogliotti2015(red))]

## define grid =============
grid <- myd[[1]]

## Attach AMC ============
amc <- vect("outputs/AMC_10m_split.shp")
zones <- rasterize(project(amc, grid),grid,field = "dz") %>% 
  as.data.table(xy=T)
zones[,`:=`(dz = factor(dz,
                        ordered = T,
                        levels = c("NorthZone_gt10", "NorthZone_lt10", 
                                   "CentralZone_gt10", "CentralZone_lt10",  
                                   "SouthZone_gt10", "SouthZone_lt10"),
                        labels = c("NorthZone_gt10", "NorthZone_lt10", 
                                   "CentralZone_gt10", "CentralZone_lt10",  
                                   "SouthZone_gt10", "SouthZone_lt10")))]

zones$dz %>% table %>% 
  as.data.table() %>% 
  knitr::kable()

d2 <- merge(d1, zones,by=c('x','y'))
d2[,`:=`(year=year(time),
         month = month(time),
         doy = yday(time))]



## Subset to dredge period =================

