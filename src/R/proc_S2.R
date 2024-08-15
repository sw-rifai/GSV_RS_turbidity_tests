pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork)
getDTthreads()

# user options ===========
product_name <- "S2-composited"
out_var <- "blue_to_green"

# region of interest ===========
amc <- vect("data/AMC/AnalysisExtentMask_GDA94z54.shp")
amc <- terra::project(amc, crs("EPSG:4326"))
amc <- amc %>% select(zone)



dir_data <- "data/gee_SAWater_turbidity/"
flist <- list.files("data/gee_SAWater_turbidity/",
  pattern = paste0("S2-composited_",out_var))
flist
flist_dates <- flist[str_detect(flist,"DATES")]
flist_r <- flist[!str_detect(flist,"DATES")]

i <- 1
tmp1 <- str_remove(flist_dates[i],"DATES_") %>% 
  str_remove(., ".csv")
tmp2 <- file.path(dir_data,flist_r[str_detect(flist_r,tmp1)]) %>% 
  sort()
fn <- function(x) crop(rast(x), amc,mask=T)
tmp3 <- lapply(tmp2, fn)
rm(tmp2); gc()


vec_dates <- fread(file.path(dir_data,flist_dates[1])) %>% 
  pull(date) %>% 
  as_date()
time(tmp3[[1]]) <- vec_dates
time(tmp3[[2]]) <- vec_dates

tmp4 <- vec_dates %>% 
  lapply(., FUN = function(x){
    terra::merge(
          tmp3[[1]][[x==time(tmp3[[1]])]],
          tmp3[[2]][[x==time(tmp3[[2]])]]
            )
  })
rm(tmp3); gc()

r <- rast(tmp4)
r1 <- r
rm(tmp4); gc(reset = T, full=T)

tmp5 <- rast(file.path(dir_data, flist_r[3]))
time(tmp5) <- fread(file.path(dir_data,flist_dates[2]))$date %>% 
  as_date()

tmp5 <- crop(tmp5, amc,mask=T)
r2 <- tmp5

r_out <- c(r1,r2)


# # system.time(
# tmp4 <- terra::merge(
#           tmp3[[1]][[x]],
#           tmp3[[2]][[x]],
#           tmp3[[3]][[x]])
# # )
# plot(tmp4,col=viridis::inferno(100))  


fp_out <- paste0(product_name,"_",out_var,"_AMC_",
         min(time(r_out)),"_",max(time(r_out)),".nc")
fp_out 
gc(reset = T, full=T)

scratch_dir <- "/home/sami/scratch"
if(dir.exists(scratch_dir)==F){
  dir.create(file.path(scratch_dir,product_name))
}
scratch_dir <- file.path(scratch_dir,product_name)

vec_time <- time(r_out)
vec_time <- sort(vec_time)
r_out[[time(r_out)==ymd("2023-11-05")]] %>% plot
fn <- function(x){
  writeCDF(r_out[[time(r_out)==x]],
    filename=file.path(scratch_dir,paste0("tmp_",x,".nc")),
    varname=out_var,
    longname=out_var,
    prec = "float",
    compression=6,
    # chunksizes = c(1069, 2904,458),
    overwrite=T, 
    verbose = F) 
}
vec_time %>% 
  lapply(., fn)

vec_time2 <- vec_time[vec_time >= ymd("2023-11-02")]

for(i in 1:length(vec_time2)){
  print(i)
  print(vec_time2[i])
  fn(vec_time2[i])
}

out_file <- file.path("data/processed_RS/",fp_out)

system(paste0("cdo -f nc4c -z zip_6 mergetime ",
  scratch_dir,
  "/*.nc ",
  out_file)
)

file.remove(list.files(scratch_dir,full.names = T))

