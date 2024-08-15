pacman::p_load(stars, sf, terra, data.table, tidyverse, lubridate)
setwd("/home/sami/srifai@gmail.com/work/research/SAWater_turbidity")
product <- "MOD09GA"
scale_factor <- 1
out_var_name <- "red"
fp <- "data/processed_RS/MOD09GA_red_AMC_2000-03-01_2023-07-01.nc"
r <- rast(fp)
vec_time <- time(r) %>% sort()

library(future.apply)
plan(multicore, workers=50)
jnk1 <- future_lapply(vec_time, FUN = function(x){
  terra::global(r[[x==time(r)]], fun="notNA") %>% 
    mutate(time = x)
})
tmp2 <- rbindlist(jnk1)
plan(sequential)

arrow::write_parquet(tmp2, sink="outputs/notNA_pix_MOD09GA_red_AMC_2000-03-01_2023-07-01.parquet")

