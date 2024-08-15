pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork)
getDTthreads() 

fn_kd490_oli <- function(bg){
  # hack factor for S2 L2A: 
  hf <- 0
  
  # https://oceancolor.gsfc.nasa.gov/resources/atbd/kd/
  vec_params <- c(-0.9054,	-1.5245,	2.2392,	-2.4777,	-1.1099)
  
  a0 <- vec_params[1]
  a1 <- vec_params[2]
  a2 <- vec_params[3]
  a3 <- vec_params[4]
  a4 <- vec_params[5]
  log10_kd490 <- a0 + a1*log10(bg) + a2*(log10(bg))**2 + a3*(log10(bg))**3 + a4*(log10(bg))**4
  kd490 <- 0.0166 + 10**log10_kd490 + hf
  return(kd490)
}
fn_kd490_oli <- Vectorize(fn_kd490_oli)

# Load 
dat <- rast("data/processed_RS/S2-composited_blue_to_green_AMC_2019-01-02_2024-03-11.nc")
dat2 <- app(dat, fun = fn_kd490_oli, cores = 50)

time(dat2) <- time(dat)
names(dat2) <- time(dat)

writeCDF(dat2, 
         filename = "data/processed_RS/S2-composited_hack-kd490_AMC_2019-01-02_2024-03-11.nc",
         varname = "kd490",
         longname = "hack kd490",
         overwrite=T)



# names(dat) <- c("x","y","time","bg")
# 
# myd <- rast("data/processed_RS/Aqua_MODIS_kd490.nc")
# time(myd)
# myd[[time(myd) == "2019-08-13"]] %>% summary # 0.06-0.9
# 
# # Add time components
# dat <- dat[is.na(bg)==F]
# dat[,`:=`(year=year(time),month=month(time))]
# vec_time <- unique(dat$time)
# 
# # calc k490
# tmp <- dat[time==ymd("2019-08-10",tz = "UTC")]
# tmp %>% 
#   filter(y > -35) %>% 
#   pull(bg) %>% `-`(0.2) %>% 
#   as.numeric() %>% 
#   fn_kd490_oli() %>% 
#   hist(100)
# ## hack factor: -0.2
# dat[,`:=`(kd490 = fn_kd490_oli(bg - 0.2))]
# 
# ref <- dat[,.(k_p99 = quantile(0.95)),by=.(x,y,month)]
# 
# 
# # curve(fn_kd490_oli,0.3,1)
# # fn_kd490_oli(0.1)
# # dat[sample(.N,1000)]$bg %>% hist(100)
# # 
# # dat$bg %>% summary
# # dat[sample(.N,10e6)]$bg %>% quantile(., 0.005)
# # fn_kd490_oli(0.43-0.1)
# # 
# # 
# # tmp %>% 
# #   filter(y > -35) %>% 
# #   ggplot(aes(x,y,fill=(fn_kd490_oli(bg-0.2))))+
# #   geom_raster()+
# #   scale_fill_viridis_c(option='H',direction=1)+
# #   coord_sf()
