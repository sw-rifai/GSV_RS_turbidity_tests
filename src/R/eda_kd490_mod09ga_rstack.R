pacman::p_load(stars, lubridate, terra, tidyverse, 
  tidyterra,
  data.table,
               mgcv, mgcViz,sf)

# Import ============
## Load raster stacks of daily MODIS kd490 (empirical)

vec_dates <- fread("data/gee_SAWater_turbidity/DATES_MOD09GA_Kd490_daily_2000-03-01_2012-12-31.csv")

ic <- rast("data/gee_SAWater_turbidity/MOD09GA_Kd490_daily_2000-03-01_2012-12-31.tif")
time(ic) <- vec_dates$date
names(ic) <- vec_dates$date 

vec_dates2 <- fread("data/gee_SAWater_turbidity/DATES_MOD09GA_Kd490_daily_2012-12-31_2023-07-01.csv")

ic2 <- rast("data/gee_SAWater_turbidity/MOD09GA_Kd490_daily_2012-12-31_2023-07-01.tif")
time(ic2) <- vec_dates2$date
names(ic2) <- vec_dates2$date

ic <- c(ic,ic2)
rm(ic2); 

r_50 <- median(ic,na.rm=T)
r_50 %>% 
  plot(col=viridis::inferno(100))
ggplot()+
  geom_spatraster(data=r_50)+
  scale_fill_viridis_c(option='B', 
    limits = c(0,1), 
    oob = scales::squish)


xy <- ic[[1]] %>% 
  as.data.table(xy=T)
xy <- xy %>% select(x,y)

dat <- terra::extract(ic, xy)
dat <- as.data.table(dat)

dat <- cbind(xy, dat)


dat <- dat %>% 
  as_tibble() %>% 
  pivot_longer(-c("x","y","ID"), 
               names_to = c("date"),
               values_to = "kd490") %>% 
  as.data.table()

bath <- vect("data/Bathymetry/DetailedModel_DEPTH.shp")
bath <- rasterize(bath, ic, field="Val_1")

amc <- vect("data/AMC/AnalysisExtentMask_GDA94z54.shp")
amc <- terra::project(amc, crs(ic))
amc <- rasterize(amc, ic, field="zone")

vv_bath <- terra::extract(bath, xy, xy=T) %>% 
  as.data.table()
vv_amc <- terra::extract(amc, xy, xy=T) %>% 
  as.data.table()
vv <- merge(vv_bath, vv_amc, by=c("x","y"))

dat <- merge(dat, vv, by=c("x","y"))
dat <- dat %>% 
  select(-starts_with("ID")) %>% 
  rename(depth = Val_1)

dat <- dat[is.na(zone)==F]
dat[,`:=`(date = ymd(date))]
dat[,`:=`(year=year(date),
  month = month(date))]

dat <- dat[,`:=`(g = case_when(
  depth < 10 ~ paste0(str_extract(zone, "\\w*"),"_lt10"),
  depth > 10 ~ paste0(str_extract(zone, "\\w*"),"_gt10")
))][is.na(g)==F]

dat[year == 2016] %>% 
  ggplot(aes(date, kd490, group=g))+
  geom_line()+
  facet_wrap(~g)

r2015 <- dat[year==2015][,.(val = median(kd490,na.rm=T)),
                           by=.(x,y,g)]
ggplot(data=r2015,aes(x,y,fill=val))+
  geom_raster()+
  scale_fill_viridis_c()+
  coord_sf()

vec_kd490_lims <- dat$kd490 %>% 
  quantile(., c(0.005,0.995),
    na.rm=T) %>% 
  unname()
dat[kd490 %between% vec_kd490_lims]$kd490 %>% hist(1000)

dat <- dat[order(g,x,y,date)] %>% 
  .[,kd490_3d := frollmean(kd490, n = 3, align = 'right',na.rm=T), 
    by = .(g,x,y)]

dat[year==2016] %>% 
  .[,.(val = median(kd490,na.rm=T)),
    by=.(date,g)] %>% 
  .[val > -50] %>% 
  .[val < 50] %>% 
  ggplot(aes(date,val,color=g))+
  geom_line()


# tmp <- (c(amc,bath)$Val_1) < 10
# bath_lt10 <- (bath$Val_1 < 10)
# bath_gt10 <- (bath$Val_1 > 10)
# as.polygons(bath_lt10) %>% plot()




ggplot()+
  geom_spatraster(data=
    bath)+
  scale_fill_viridis_c()


gc()

ic_u <- mean(ic,na.rm=T)
ic_50 <- median(ic,na.rm=T)
ic_sd <- stdev(ic,na.rm=T)

ggplot()+
  geom_spatraster(data=ic_u)+
  scale_fill_viridis_c(limits=c(-2,1))
ggplot()+
  geom_spatraster(data=ic_50)+
  scale_fill_viridis_c(limits=c(-0.5,0),oob=scales::squish)
ggplot()+
  geom_spatraster(data=ic_sd)+
  scale_fill_viridis_c(
    # limits=c(0,10)
    )

fn <- function(r){
  out <- as.data.table(r,xy=T)
  out$date <- names(r)
  return(out)
}

ic[[1:2]]

dat <- lapply(ic[[1:3]], fn)

dat <- stars::st_as_stars(ic) %>% 
   as.data.table(xy=T)

names(dat) <- c("x","y","")