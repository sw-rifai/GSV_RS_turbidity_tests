pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork)

# OPTIONS ===================================
var_name <- "red" 
sensor2 <- "MCD43A4-qaMasked"

fp <- "data/processed_RS/MCD43A4-qaMasked_red_AMC_2000-02-24_2023-10-22.nc"
# ===========================================


fn_dogliotti <- function(p_w){
  # Dogliotti et al 2015, Remote Sensing of Environment
  # Notes: See ACOLITE for most recent parameter calibrations
  A_t <- 378.46
  C <- 0.19905
  out <- A_t * (p_w) / 
    (1 - ((p_w)/ C))
  return(out)
}
fn_dogliotti <- Vectorize(fn_dogliotti)

r <- rast(fp)

tmp <- list()
for(i in 1:100){
  tmp[[i]] <- fn_dogliotti(r[[i]])
}

tmp2 <- lapply(tmp, rast)

tmp3 <- rast(tmp2)
tmp4 <- mean(tmp3,na.rm=T)

plot(tmp3[[50]], col=inferno(100))

tmp <- fn_dogliotti(r[[1]])
rast(tmp) %>% plot
plot(tmp)

r[[1000]] %>% is.null
r[[1]] %>% is.null
r[[1]] %>% is.null

tmp <- r[[1:100]]


fn_dogliotti(runif(100,0,0.1))
fn_dogliotti(NULL)

d <- lapp(r[[1:100]], fn)
