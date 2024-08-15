pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork)

scale_factor <- 1 # 0.0001;

# Import ============ 

# Torrens river discharge ------------------------
flow <- fread("data/daily_discharge_torrens_river_seaview_road_bridge.csv")
names(flow) <- c("date","time","dv")
flow[,`:=`(date = dmy(date))]
flow[,`:=`(year=year(date),month=month(date))]
flow <- flow[order(date)][,`:=`(dv_max5d = frollapply(dv,5,max,na.rm=T,align='right'))]
flow <- flow[order(date)][,`:=`(dv_max15d = frollapply(dv,15,max,na.rm=T,align='right'))]
flow <- flow[order(date)][,`:=`(dv_max30d = frollapply(dv,30,max,na.rm=T,align='right'))]


## Load and collate raster stacks of daily MODIS SMI
vec_dates <- fread("data/gee_SAWater_turbidity/DATES_MOD09GA_kd490_daily_2000-03-01_2012-12-31.csv")
vec_dates2 <- fread("data/gee_SAWater_turbidity/DATES_MOD09GA_kd490_daily_2012-12-31_2023-07-01.csv")

ic <- rast("data/gee_SAWater_turbidity/MOD09GA_kd490_daily_2000-03-01_2012-12-31.tif")
ic2 <- rast("data/gee_SAWater_turbidity/MOD09GA_kd490_daily_2012-12-31_2023-07-01.tif")

time(ic) <- vec_dates$date
time(ic2) <- vec_dates2$date
names(ic) <- vec_dates$date
names(ic2) <- vec_dates2$date

ic <- c(ic,ic2)
rm(ic2)

## examine LT estimates
ic_u <- mean(ic,na.rm=T)
ic_50 <- median(ic,na.rm=T)
ic_sd <- terra::stdev(ic,na.rm=T)

plot(ic_u, col = inferno(10))
plot(ic_50, col = inferno(10))
plot(ic_sd, col = inferno(10))

## Define coords
xy <- ic_u %>% 
  as.data.table(xy=T) %>% 
  select(x,y)

## extract depth and AMC zones
bath <- vect("data/Bathymetry/DetailedModel_DEPTH.shp")
bath <- rasterize(bath, ic, field="Val_1")
names(bath) <- "depth"

amc <- vect("data/AMC/AnalysisExtentMask_GDA94z54.shp")
amc <- terra::project(amc, crs(ic))
amc <- rasterize(amc, ic, field="zone")

vv_bath <- terra::extract(bath, xy[,.(x,y)], xy=T) %>% 
  as.data.table()
vv_amc <- terra::extract(amc, xy[,.(x,y)], xy=T) %>% 
  as.data.table()

tmp <- terra::extract(ic, xy)
tmp <- cbind(xy,tmp)

dat <- cbind(tmp, vv_bath %>% select(depth),
             vv_amc %>% select(zone))


dat <- dat %>% 
  as_tibble() %>% 
  pivot_longer(-c("x","y","ID","depth","zone"), 
               names_to = c("date"),
               values_to = "kd490") %>% 
  as.data.table()

##

dat <- dat[is.na(zone)==F]
dat[,`:=`(date = ymd(date))]
dat[,`:=`(year=year(date),
          month = month(date))]

dat <- dat[,`:=`(g = case_when(
  depth < 10 ~ paste0(str_extract(zone, "\\w*"),"_lt10"),
  depth > 10 ~ paste0(str_extract(zone, "\\w*"),"_gt10")
))][is.na(g)==F]

gc()

unique(dat[,.(g,x,y)]) %>% pull(g) %>% table()

dat <- dat[order(g,x,y,date)]
# dat <- dat[,`:=`(kd490_max5d = frollapply(kd490, n=5, 
#                                max, fill=NA, align='right',na.rm=T)),
#    by=.(g)]





# ss <- dat[year==2016][month %between% c(8,12)][kd490_max5d > 0][is.na(g)==F] %>%
#   .[,.(kd490_max5d = median(kd490_max5d,na.rm=T), 
#        nobs = sum(is.na(kd490_max5d)==F)),
#     by=.(date,g)]
# ss <- ss[,`:=`(year=year(date),
#          month= month(date))]

sm <- dat[,.(kd49050 = median(kd490,na.rm=T), 
             nobs = sum(is.na(kd490)==F)),
          by=.(g,date,year,month)]
sm <- sm[order(g,date)]
sm <- sm[,`:=`(kd490_max5d = frollapply(kd49050, n=5, 
                         max, fill=NA, align='right',na.rm=T)),
           by=.(g)]

# ss <- dat[year==2016][month %between% c(8,12)][is.na(g)==F] %>%
#   .[,.(kd490_max5d = median(kd490_max5d,na.rm=T),
#        nobs = sum(is.na(kd490_max5d)==F)),
#     by=.(date,g)]


merge(flow[year==2016], sm, 
      by=c("date","month","year")) %>% 
  ggplot(aes(dv_max5d, kd490_max5d, color = g))+
  geom_point()+
  geom_smooth(se = F, 
              method = 'lm')


(p1 <- flow[year==2016][month %between% c(8,12)] %>% 
  ggplot(aes(date, dv_max5d))+
  # geom_point()+ 
  geom_line() + 
  labs(x = NULL, 
       y = "Rolling 5-day Max Daily Volume (ML)")+ 
  scale_x_date(date_labels = "%Y-%m-%d") + 
  theme_linedraw() + 
  theme(panel.grid = element_blank()))



(p2 <- sm[year==2016][month %between% c(8,12)] %>% 
  ggplot(aes(date, kd490_max5d*scale_factor,color=g))+
  geom_line()+
  geom_point(aes(size = nobs), 
             alpha = 0.85, 
             shape = 21) + 
  scale_color_discrete_c4a_cat(palette = 'tol.muted') +
  labs( y = "Rolling 5-day max of zonal kd490", 
        x = NULL, 
        size = "Number of pixel obs.",
        color = "Zone and depth") + 
  scale_x_date(date_labels = "%Y-%m-%d") +
  theme_linedraw() + 
  theme(panel.grid = element_blank(), 
        legend.position = 'bottom'))

(p_out <- p1/p2)

ggsave(p_out,
       filename = "figures/figure_eda_mod09ga_kd490-5dmax_2016.png",
       width = 30,
       height = 20,
       units = 'cm',
       dpi = 350)

# 
# unique(dat[,.(x,y,g)]) %>% 
#   pull(g) %>% 
#   table()
# 
# 
# flow[dv==max(dv)]
# 
# dat[year==2016] %>% 
#   ggplot(aes(date, poc,color=g))+
#   geom_point()
#   geom_smooth(se=F, 
#               method = 'gam',
#               formula = y~s(x,bs='ad'))
# 
# 
# dat[year==2016] %>% 
#   ggplot(aes(date, poc,color=g))+
#   geom_line()
