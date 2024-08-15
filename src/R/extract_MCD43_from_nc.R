pacman::p_load(stars, sf, terra, data.table, tidyverse, lubridate)
library(exactextractr)

setwd("/home/sami/srifai@gmail.com/work/research/SAWater_turbidity")

product <- "MCD43A4"
scale_factor <- 1
out_var_name <- "blue_to_green"
fp <- "data/processed_RS/MCD43A4-qaMasked_blue_to_green_AMC_2000-02-24_2023-10-22.nc"
r <- rast(fp)

# plot(r[[1]],col=viridis(100))

# dir_data <- "/media/sami/drive_a/landsat_marine/AMC/"
# flist <- list.files(dir_data,
#                     recursive = T, 
#                     pattern = "AR_BAND4")
# flist

amc <- sf::read_sf('outputs/AMC_10m_split.gpkg')
amc <- vect(amc)
crs_r <- r %>% crs()
amc <- project(amc, crs_r)
amc_s <- st_as_sf(amc)

vec_time <- time(r) %>% sort()

t1 <- arrow::read_parquet("outputs/notNA_pix_MCD43A4-qaMasked_red_AMC_2000-02-24_2023-10-22.parquet")
t1 <- t1[notNA >= 250]
vec_time <- vec_time[vec_time %in% t1$time]

fn_get_vals <- function(in_time){
  # tmp_time <- str_extract(fp,"_\\d{8}_") %>% 
  #   str_remove(., "_") %>% 
  #   str_remove(.,"_") %>% 
  #   ymd()
  
  tmp1 <- r[[time(r)==in_time]]
  tmp1 <- crop(tmp1,amc,mask=T)
  
  tmp2 <- exactextractr::exact_extract(tmp1, amc_s, fun = c('mean', 'median', 'count','stdev','quantile'), 
                        quantiles=c(0.01, 0.05, 0.25,0.75, 0.95,0.99),
                        append_cols = 'dz') %>% setDT()
  tmp3 <- cbind(tmp2$dz,tmp2[,2:ncol(tmp2)]*scale_factor)
  names(tmp3) <- names(tmp2)
  tmp3$time <- in_time
  # tmp3$src_file <- fp
  return(tmp3)
}

library(future.apply)
plan(multicore, workers=50)
out1 <- future_lapply(vec_time, FUN = fn_get_vals)
out2 <- rbindlist(out1)
plan(sequential)



fp_out <- paste0(product,"_AMC_",out_var_name,".parquet")
fp_out
arrow::write_parquet(out2, 
                     sink=file.path("outputs",fp_out),
                     compression = 'snappy')


# jnk1 <- arrow::read_parquet("outputs/MCD43A4_AMC_red.parquet")

# library(future.apply)
# tmp <- vec_time[1:10] %>% 
#   lapply(., 
#          function(x){
#            global(r[[time(r)==x]], fun="notNA")
#          }
#   )
# 
# tmp
# 
# plan(sequential)
# s_list <- vec_time[1:10]
# jnk1 <- future_lapply(s_list, FUN = function(x){
#     terra::global(r[[x==time(r)]], fun="notNA") %>% 
#     mutate(time = x)
#   })
# jnk1
# 
# dat
# # dat2 <- dat %>% 
# #   rbindlist()

# jnk1 %>% 
#   ggplot(aes(time, median,color=dz))+
#   geom_line()+
#   facet_wrap(~dz,nrow=3)

