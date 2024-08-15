pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork)

# Import ============

# Torrens river discharge ------------------------
flow <- fread("data/daily_discharge_torrens_river_seaview_road_bridge.csv")
names(flow) <- c("date","time","dv")
flow[,`:=`(date = dmy(date))]
flow[,`:=`(year=year(date),month=month(date))]


## Load and collate raster stacks of daily MODIS SMI
vec_dates <- fread("data/gee_SAWater_turbidity/DATES_L3SMI_poc_daily_2002-07-03_2012-12-31.csv")
vec_dates2 <- fread("data/gee_SAWater_turbidity/DATES_L3SMI_poc_daily_2012-12-31_2023-07-01.csv")

ic <- rast("data/gee_SAWater_turbidity/L3SMI_poc_daily_2002-07-03_2012-12-31.tif")
ic2 <- rast("data/gee_SAWater_turbidity/L3SMI_poc_daily_2012-12-31_2023-07-01.tif")

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

tmp <- terra::extract(ic, xy)
tmp <- cbind(xy,tmp) %>% 
  as_tibble() %>% 
  pivot_longer(-c("x","y","ID"), 
               names_to = c("date"),
               values_to = "poc") %>% 
  as.data.table()


## extract depth and AMC zones
bath <- vect("data/Bathymetry/DetailedModel_DEPTH.shp")
bath <- rasterize(bath, ic, field="Val_1")
names(bath) <- "depth"

amc <- vect("data/AMC/AnalysisExtentMask_GDA94z54.shp")
amc <- terra::project(amc, crs(ic))
amc <- rasterize(amc, ic, field="zone")

vv_bath <- terra::extract(bath, tmp[,.(x,y)], xy=T) %>% 
  as.data.table()
vv_amc <- terra::extract(amc, tmp[,.(x,y)], xy=T) %>% 
  as.data.table()

dat <- cbind(tmp, vv_bath %>% select(depth),
             vv_amc %>% select(zone))

dat <- dat[is.na(zone)==F]
dat[,`:=`(date = ymd(date))]
dat[,`:=`(year=year(date),
          month = month(date))]

dat <- dat[,`:=`(g = case_when(
  depth < 10 ~ paste0(str_extract(zone, "\\w*"),"_lt10"),
  depth > 10 ~ paste0(str_extract(zone, "\\w*"),"_gt10")
))][is.na(g)==F]


sm <- dat[order(g,x,y,date)][,`:=`(poc_max5d = frollapply(poc, n=5, max, fill=NA, align='right',na.rm=T)),
   by=.(g)]


p1 <- flow[year==2016][month %between% c(8,12)] %>% 
  ggplot(aes(date, dv))+
  # geom_point()+ 
  geom_line() + 
  labs(x = NULL, 
       y = "Daily Volume (ML)")+ 
  scale_x_date(date_labels = "%Y-%m-%d") + 
  theme_linedraw() + 
  theme(panel.grid = element_blank())



p2 <- sm[year==2016][month %between% c(8,12)][poc_max5d > 0][is.na(g)==F] %>%
  .[,.(poc = median(poc,na.rm=T)),
    by=.(date,g)] %>% 
  ggplot(aes(date, poc,color=g))+
  geom_line()+
  geom_point(size = 2) + 
  scale_color_discrete_c4a_cat(palette = 'tol.muted') +
  labs( y = "Daily POC", 
        x = NULL, 
        color = "Zone and depth") + 
  scale_x_date(date_labels = "%Y-%m-%d") +
  theme_linedraw() + 
  theme(panel.grid = element_blank())

p1/p2

ggsave("figures/figure_eda_l3smi_poc.png",
       width = 25,
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
