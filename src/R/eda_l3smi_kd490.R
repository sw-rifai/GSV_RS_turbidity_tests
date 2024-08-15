pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork) 

scale_factor <- 0.0001;

# Import ============ 

# Torrens river discharge ------------------------
flow <- fread("data/daily_discharge_torrens_river_seaview_road_bridge.csv")
names(flow) <- c("date","time","dv")
flow[,`:=`(date = dmy(date))]
flow[,`:=`(year=year(date),month=month(date))]
flow <- flow[order(date)][,`:=`(dv_max5d = frollapply(dv,5,max,na.rm=T,align='right'))]
flow <- flow[order(date)][,`:=`(dv_max15d = frollapply(dv,15,max,na.rm=T,align='right'))]
flow <- flow[order(date)][,`:=`(dv_max30d = frollapply(dv,30,max,na.rm=T,align='right'))]


## Load and collate raster stacks of daily MODIS NIR
vec_dates <- fread("data/gee_SAWater_turbidity/DATES_L3SMI_combined_kd490_mod_daily_2002-07-03_2012-12-31.csv")
vec_dates2 <- fread("data/gee_SAWater_turbidity/DATES_L3SMI_combined_kd490_mod_daily_2012-12-31_2023-07-01.csv")

ic <- rast("data/gee_SAWater_turbidity/L3SMI_combined_kd490_mod_daily_2002-07-03_2012-12-31.tif")
ic2 <- rast("data/gee_SAWater_turbidity/L3SMI_combined_kd490_mod_daily_2012-12-31_2023-07-01.tif")

time(ic) <- vec_dates$date
time(ic2) <- vec_dates2$date
names(ic) <- vec_dates$date
names(ic2) <- vec_dates2$date

ic <- c(ic,ic2)
rm(ic2)



## examine LT estimates
ic_u <- mean(ic,na.rm=T)
# ic_50 <- median(ic,na.rm=T)
# ic_sd <- terra::stdev(ic,na.rm=T)

# plot(ic_u, col = inferno(10))
# plot(ic_50, col = inferno(10))
# plot(ic_sd, col = inferno(10))

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

tmp_red <- terra::extract(ic, xy)
tmp_red <- cbind(xy,tmp_red)

dat <- cbind(tmp_red, vv_bath %>% select(depth),
             vv_amc %>% select(zone))


cdat <- dat %>% 
  as_tibble() %>% 
  pivot_longer(-c("x","y","ID","depth","zone"), 
               names_to = 'date',
               values_to = "kd490") %>% 
  as.data.table()

cdat <- cdat[is.na(zone)==F]
cdat[,`:=`(date = lubridate::ymd(date))]


# library(stars)
# cdat <- c(st_as_stars(nir) %>% 
#   set_names("nir"),
#   st_as_stars(red) %>% 
#   set_names("red")) %>% 
#   as.data.table(xy=T)
# 
# cdat <- merge(cdat, 
#       merge(vv_bath,vv_amc,by=c("x","y","ID")),
#       by=c("x","y"))


# cdat[,`:=`(ndvi = (nir-red)/(nir+red))]
# cdat[,`:=`(r_nr = red/(nir+red))]
# cdat[,`:=`(n_nr = nir/(nir+red))]


cdat <- cdat[,`:=`(g = case_when(
  depth < 10 ~ paste0(str_extract(zone, "\\w*"),"_lt10"),
  depth > 10 ~ paste0(str_extract(zone, "\\w*"),"_gt10")
))][is.na(g)==F]


# cdat <- cdat %>% rename(date = time)
# cdat[,`:=`(date = date(date))]

sm <- cdat[,.(
              # ndvi50 = median(ndvi,na.rm=T), 
              # r_nr50 = median(r_nr,na.rm=T),
              # n_nr50 = median(n_nr,na.rm=T), 
              val95  = quantile(kd490,0.95, na.rm=T),
              val50  = median(kd490,na.rm=T),
              val_sd = sd(kd490,na.rm=T)
            # n50 = median(nir,na.rm=T)
            ), 
           by=.(date,g)]


sm[,`:=`(date = ymd(date))]
sm[,`:=`(year=year(date),
          month = month(date))]

sm[year==2016][month %in% c(8:12)] %>% 
  ggplot(aes(date, val_sd, color = g)) +
  geom_line()
# cdat[sample(1e5)] %>% ggplot(aes(nir,red))+geom_point()+geom_abline(col='red')




sm <- cdat[order(g,x,y,date)][,
                 `:=`(kd_max5d = frollapply(kd490, n=5, max, fill=NA, align='right',na.rm=T)),
                   by=.(g)]
sm[,`:=`(year=year(date),
         month= month(date))]

ss <- sm %>% 
  .[year==2016] %>%
  .[month %between% c(8,12)] %>%
  .[kd_max5d > 0] %>%
  .[is.na(g)==F] %>%
  .[,.(kd_max5d = median(kd_max5d,na.rm=T), 
       kd490 = median(kd490,na.rm=T), 
       nobs = sum(is.na(kd_max5d)==F)),
    by=.(date,g)]
ss[,`:=`(date = date(date))]
ss <- ss[,`:=`(year=year(date),
               month= month(date))]


# # Plotting ========================================================
# 
# ## AMC zone plot ===============
# oz_poly <- rnaturalearth::ne_countries(country = "Australia",
#                                        returnclass = "sf", 
#                                        scale = 'large')
# amc <- sf::read_sf("data/AMC/AnalysisExtentMask_GDA94z54.shp")
# amc <- sf::st_transform(amc, sf::st_crs(oz_poly))
# amc <- amc %>% mutate(
#   fzone = factor(zone, ordered = T,
#                  levels = c("North Zone","Central Zone","South Zone"),
#                  labels = c("North Zone","Central Zone","South Zone")))
# 
# vec_lims <- st_bbox(st_buffer(amc,10000))
# 
# ggplot()+
#   geom_sf(data = oz_poly,
#           fill = "grey30")+ 
# geom_sf(data = amc, aes(fill = fzone)) +
#   coord_sf(xlim = vec_lims[c("xmin","xmax")],
#            ylim = vec_lims[c("ymin","ymax")], 
#            crs = st_crs(amc), 
#            ndiscr = 3,
#            lims_method = "box") + 
#   scale_x_continuous(breaks = seq(
#     from = round(vec_lims["xmin"],digits=1), 
#     to = round(vec_lims['xmax'],digits=1), 
#     by = 0.3
#   )) +
#   scale_fill_viridis_d(option = "C", end = 0.9) + 
#   labs(fill = "Zone") + 
#   theme(panel.background = element_rect(fill='lightblue'), 
#         panel.grid = element_blank(), 
#         legend.position = 'bottom')
# ggsave(filename = "figures/AMC_zones.png",
#        width = 11,
#        height = 15,
#        units = 'cm',
#        dpi = 350)


## Map of Retrievals during October
cdat[date %in% seq(ymd("2016-10-01"),ymd("2016-11-01"),by="1 day")] %>% 
  ggplot(aes(x,y,fill=kd490))+
  geom_raster()+
  coord_sf()+
  scale_fill_viridis_c(option = 'B', 
                       begin = 0.1,
                       na.value = "grey70",
                       # limits = c(0,0.2),
                       oob = scales::squish) + 
  labs(fill = "kd490", 
       title = "L3SMI Adelaide Metro. Coastal Zones, October 2016", 
       x = NULL,
       y = NULL) +
  facet_wrap(~date, 
             nrow = 4) + 
  theme_linedraw()+
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(size = 7),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = 'grey20'), 
        legend.key.width = unit(1, 'cm'),
        legend.position = 'bottom')
ggsave(filename = "figures/L3SMI_2016-10_map_kd490.png",
       width = 17,
       height = 25,
       units = 'cm',
       dpi = 350)


ss %>% 
  ggplot(aes(date,kd490,color=g))+
  geom_point()


merge(flow, ss, 
      by=c("date","month","year")) %>% 
  ggplot(aes(dv_max5d, 
             kd_max5d, color = g))+
  geom_point()+
  geom_smooth(se = F, 
              method = 'lm') + 
  facet_wrap(~g)

merge(flow[year==2016], ss, 
      by=c("date","month","year")) %>% 
  ggplot(aes(dv_max5d, red, color = g))+
  # geom_point()+
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



(p2 <- ss %>% 
    ggplot(aes(date, kd_max5d*scale_factor,color=g))+
    geom_line()+
    geom_point(aes(size = nobs), 
               alpha = 0.85, 
               shape = 21) + 
    scale_color_discrete_c4a_cat(palette = 'tol.muted') +
    labs( y = "5-day Max Red Reflectance", 
          x = NULL, 
          size = "Number of pixel obs.",
          color = "Zone and depth") + 
    scale_x_date(date_labels = "%Y-%m-%d") +
    theme_linedraw() + 
    theme(panel.grid = element_blank(), 
          legend.position = 'bottom'))

(p_out <- p1/p2)

ggsave(p_out,
       filename = "figures/figure_eda_mcd43a4_red.png",
       width = 30,
       height = 20,
       units = 'cm',
       dpi = 350)



## 
ss2 <- sm %>% 
  # .[year==2016] %>%
  # .[month %between% c(8,12)] %>%
  .[kd_max5d > 0] %>%
  .[is.na(g)==F] %>%
  .[,.(kd_max5d = median(kd_max5d,na.rm=T), 
       kd490 = median(kd490,na.rm=T), 
       nobs = sum(is.na(kd_max5d)==F)),
    by=.(date,g)]
ss2[,`:=`(date = date(date))]
ss2[,`:=`(ddate = decimal_date(date))]
ss2 <- ss2[,`:=`(year=year(date),
               month= month(date))]



m1 <- bam(kd_max5d ~ 
            te(x,y,ddate,k=c(3,3,30)), 
          data = sm[kd_max5d>0] %>% 
            mutate(ddate = decimal_date(date)), 
          family=Gamma(link='log'),
          select = T,
          discrete = F)
summary(m1)
plot(m1)
getViz(m1) %>% plot() %>% 
  print(pages = 1)

pdat <- expand_grid(unique(sm[,.(x,y)]), 
            ddate = unique(sm %>% mutate(ddate=decimal_date(date)) %>% pull(ddate))) %>% 
  as.data.table()
pdat <- pdat %>% 
  mutate(pred = predict(m1,newdata=.,type='response'))

sm[,`:=`(ddate = decimal_date(date))]

pdat[,.(val = mean(pred)),by=ddate] %>% 
  ggplot(aes(ddate,val))+
  geom_line(data=sm[,.(val = mean(kd490,na.rm=T)),by=ddate],
            aes(ddate,val),col='red') + 
  geom_line()
  



sm %>% 
  .[year==2016] %>%
  .[month %between% c(8,12)] %>%
  .[kd_max5d > 0] %>%
  .[is.na(g)==F] %>%
  .[,.(kd_max5d = median(kd_max5d,na.rm=T), 
       kd490 = median(kd490,na.rm=T), 
       nobs = sum(is.na(kd_max5d)==F)),
    by=.(date,g)]



