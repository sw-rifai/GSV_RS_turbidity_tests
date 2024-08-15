pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork)

scale_factor <- 0.0001;
peak_flow_date <- flow[dv==max(dv)]$date

# Import ============ 

# Torrens river discharge ------------------------
flow <- fread("data/daily_discharge_torrens_river_seaview_road_bridge.csv")
names(flow) <- c("date","time","dv")
flow[,`:=`(date = dmy(date))]
flow[,`:=`(year=year(date),month=month(date))]
flow <- flow[order(date)][,`:=`(dv_max5d = frollapply(dv,5,max,na.rm=T,align='right'))]
flow <- flow[order(date)][,`:=`(dv_max15d = frollapply(dv,15,max,na.rm=T,align='right'))]
flow <- flow[order(date)][,`:=`(dv_max30d = frollapply(dv,30,max,na.rm=T,align='right'))]


## Load and collate raster stacks of daily MODIS NIR -------------------------
vec_dates <- fread("data/gee_SAWater_turbidity/DATES_MCD43A4_Nadir_Reflectance_Band2_daily_2002-02-04_2012-12-31_c2e7c3453c422a1824742911bb6864ed.csv")
vec_dates2 <- fread("data/gee_SAWater_turbidity/DATES_MCD43A4_Nadir_Reflectance_Band2_daily_2012-12-31_2023-07-11_c2e7c3453c422a1824742911bb6864ed.csv")

ic <- rast("data/gee_SAWater_turbidity/MCD43A4_Nadir_Reflectance_Band2_daily_2002-02-04_2012-12-31_c2e7c3453c422a1824742911bb6864ed.tif")
ic2 <- rast("data/gee_SAWater_turbidity/MCD43A4_Nadir_Reflectance_Band2_daily_2012-12-31_2023-07-11_c2e7c3453c422a1824742911bb6864ed.tif")

time(ic) <- vec_dates$date
time(ic2) <- vec_dates2$date
names(ic) <- vec_dates$date
names(ic2) <- vec_dates2$date

ic <- c(ic,ic2)
rm(ic2)

## examine LT estimates
ic_u <- mean(ic,na.rm=T)

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

tmp_nir <- terra::extract(ic, xy)
tmp_nir <- cbind(xy,tmp_nir)

dat_nir <- cbind(tmp_nir, vv_bath %>% select(depth),
             vv_amc %>% select(zone))


cdat_nir <- dat_nir %>% 
  as_tibble() %>% 
  pivot_longer(-c("x","y","ID","depth","zone"), 
               names_to = 'date',
               values_to = "nir") %>% 
  as.data.table()
## END SECTION ------------------------


## Load and collate raster stacks of daily MODIS RED -------------------------
vec_dates <- fread("data/gee_SAWater_turbidity/DATES_MCD43A4_Nadir_Reflectance_Band1_daily_2002-02-04_2012-12-31_61e3bb679d81da01099535c47701ca79.csv")
vec_dates2 <- fread("data/gee_SAWater_turbidity/DATES_MCD43A4_Nadir_Reflectance_Band1_daily_2012-12-31_2023-07-11_61e3bb679d81da01099535c47701ca79.csv")

ic <- rast("data/gee_SAWater_turbidity/MCD43A4_Nadir_Reflectance_Band1_daily_2002-02-04_2012-12-31_61e3bb679d81da01099535c47701ca79.tif")
ic2 <- rast("data/gee_SAWater_turbidity/MCD43A4_Nadir_Reflectance_Band1_daily_2012-12-31_2023-07-11_61e3bb679d81da01099535c47701ca79.tif")

time(ic) <- vec_dates$date
time(ic2) <- vec_dates2$date
names(ic) <- vec_dates$date
names(ic2) <- vec_dates2$date

ic <- c(ic,ic2)
rm(ic2)

ic_u <- mean(ic,na.rm=T)

## Define coords
xy <- ic_u %>% 
  as.data.table(xy=T) %>% 
  select(x,y)

tmp_red <- terra::extract(ic, xy)
tmp_red <- cbind(xy,tmp_red)

# dat_red <- cbind(tmp_red, vv_bath %>% select(depth),
#              vv_amc %>% select(zone))


cdat_red <- tmp_red %>% 
  as_tibble() %>% 
  pivot_longer(-c("x","y","ID"), 
               names_to = 'date',
               values_to = "red") %>% 
  as.data.table()
## END SECTION -------------

## merge stuff
system.time(
cdat <- merge(cdat_red, cdat_nir, by=c("x","y","ID","date"))
)
cdat <- cdat[is.na(depth)==F][is.na(zone)==F]
cdat[,`:=`(date = lubridate::ymd(date))]

## groups
cdat <- cdat[,`:=`(g = case_when(
  depth < 5 ~ paste0(str_extract(zone, "\\w*"),"_lt5"),
  depth < 10 ~ paste0(str_extract(zone, "\\w*"),"_lt10"),
  depth > 10 ~ paste0(str_extract(zone, "\\w*"),"_gt10")
))][is.na(g)==F]


## calc metrics
cdat[,`:=`(nir = nir*0.0001, 
           red = red*0.0001)]
cdat[,`:=`(nir_red = nir/red)]
cdat[,`:=`(tss = 1-nir_red)]
cdat[,`:=`(ndvi = (nir-red)/(nir+red))]


cdat[sample(.N, 1e6)] %>% select(red,nir) %>% drop_na() %>% 
  cor()
cdat[sample(.N, 1e6)][red > 0] %>% pull(tss) %>% 
  quantile(., c(0.01,0.99), na.rm=T)
  hist(100)


cdat <- cdat[order(g,x,y,date)][,
          `:=`(nir_max5d = frollapply(nir, n=5, max, fill=NA, align='right',na.rm=T),
               red_max5d = frollapply(red, n=5, max, fill=NA, align='right',na.rm=T)),
          by=.(g)]
  

# cdat <- cdat %>% rename(date = time)
# cdat[,`:=`(date = date(date))]

sm <- cdat[,.(
  nobs = sum(is.na(red)==F), 
  red_max5d95 = quantile(red_max5d, 0.95,na.rm=T),
  red_max5d50 = median(red_max5d,na.rm=T), 
  
  tss50  = median(tss,na.rm=T),
  tss95  = quantile(tss,0.95, na.rm=T),
  tss_sd = sd(tss,na.rm=T),
  red50  = median(red,na.rm=T),
  red95  = quantile(red,0.95, na.rm=T),
  red_sd = sd(red,na.rm=T), 
  nir50  = median(nir,na.rm=T),
  nir95  = quantile(nir,0.95, na.rm=T),
  nir_sd = sd(nir,na.rm=T)
  # n50 = median(nir,na.rm=T)
), 
by=.(date,g)][,`:=`(date = ymd(date))][,`:=`(year=year(date),
                                             month = month(date))]

sm_ref <- sm[,.(red50_u = mean(red50,na.rm=T), 
                red95_u = mean(red95,na.rm=T),
                nir95_u = mean(nir95,na.rm=T), 
                nir50_u = mean(nir50,na.rm=T),
                tss95_u = mean(tss95,na.rm=T)), 
             by=.(month,g)]



(p1 <- flow[year==2016][month %between% c(8,12)] %>% 
    ggplot(aes(date, dv_max5d))+
    # geom_point()+ 
    geom_line() + 
    labs(x = NULL, 
         y = "Rolling 5-day Max Daily Volume (ML)")+ 
    scale_x_date(date_labels = "%Y-%m-%d") + 
    theme_linedraw() + 
    theme(panel.grid = element_blank()))

(p2 <- sm[str_detect(g, "lt")][year==2016][month %between% c(8,12)] %>% 
    select(date,g, red_max5d95, red_max5d50) %>% 
    pivot_longer(-c("date","g")) %>% 
    ggplot(aes(date, value,color=g))+
    geom_vline(aes(xintercept = peak_flow_date), 
               lty = 3, col = 'grey30') + 
    geom_line()+
    # geom_point(aes(size = nobs), 
    #            alpha = 0.85, 
    #            shape = 21) + 
    scale_color_discrete_c4a_cat(palette = 'tol.muted') +
    labs( y = "5-day Max nir Reflectance", 
          x = NULL, 
          size = "Number of pixel obs.",
          color = "Zone and depth") + 
    scale_x_date(date_labels = "%Y-%m-%d") +
    facet_wrap(~name, ncol = 1, scales = 'free') + 
    theme_linedraw() + 
    theme(panel.grid = element_blank(), 
          legend.position = 'bottom'))

(p_out <- p1/p2 + plot_layout(heights = c(1,2.5)))

ggsave(p_out,
       filename = "figures/figure_eda_mcd43a4_red5d95_red5d50.png",
       width = 30,
       height = 20,
       units = 'cm',
       dpi = 350)


mdat <- cdat %>% 
  # [date %between% c(ymd("2016-10-01"),ymd("2017-10-01"))] %>% 
  mutate(month = month(date), 
         year = year(date)) %>% 
  .[,.(red = quantile(red, 0.95, na.rm=T), 
       nobs = sum(is.na(red)==F)), 
    by=.(year,month,x,y,depth,zone)] %>% 
  mutate(date = ymd(paste(year,month,1))) 
mdat_ref <- mdat[,.(red_u = mean(red, na.rm=T)), 
                 by=.(zone, depth, month)]

merge(mdat[date %between% c(ymd("2016-08-01"),ymd("2017-07-01"))], 
      mdat_ref, 
      by=c("zone","month","depth")) %>% 
  mutate(red_anom = red - red_u) %>% 
  .[depth < 20] %>% 
  ggplot(aes(depth, red_anom, color=zone))+
  geom_hline(aes(yintercept = 0)) + 
  geom_smooth(aes(weight = nobs),
              se = F, 
              method = 'bam',
              formula = y~s(x, bs='ts', k=10), 
              method.args = list(select=T, 
                                 discrete = T)) + 
  scale_color_discrete_c4a_cat(palette = 'tol.muted') +
  labs(x = "Depth (m)",
       y= "Monthly Anomaly of 95th percentile of Red Reflectance", 
       color = "Zone") + 
  facet_wrap(~date) + 
  theme_linedraw()+
  theme(panel.grid = element_blank(), 
        legend.position = 'bottom',
        strip.background = element_blank(),
        strip.text = element_text(color='black'))

ggsave(
       filename = "figures/figure_eda_mcd43a4_red95p_anom.png",
       width = 20,
       height = 15,
       units = 'cm',
       dpi = 350)


# 
# 
# sm[date %between% c(ymd("2016-08-01"),ymd("2017-01-01"))][str_detect(g,'_lt5')] %>% 
#   merge(., sm_ref, by=c("month","g")) %>% 
#   ggplot(aes(date,
#              100*(tss95 - tss95_u)/tss95_u,
#              # red95, 
#              # 100*(red_max5d95-red95_u)/red95_u,
#              color = g)) +
#   geom_line() + 
#   scale_color_discrete_c4a_cat() + 
#   theme_linedraw() +
#   theme()
# 
# sm[date %between% c(ymd("2015-08-01"),ymd("2016-08-01"))][str_detect(g,'_lt5')] %>% 
#   ggplot(aes(date, red95, color = g)) +
#   geom_line() + 
#   scale_color_discrete_c4a_cat() + 
#   theme_linedraw() +
#   theme()
# 
# # cdat[sample(1e5)] %>% ggplot(aes(nir,nir))+geom_point()+geom_abline(col='nir')
# 
# 
# 
# 
# sm[,`:=`(year=year(date),
#          month= month(date))]
# 
# ss <- sm %>% 
#   .[year==2016] %>%
#   .[month %between% c(8,12)] %>%
#   .[nir_max5d > 0] %>%
#   .[is.na(g)==F] %>%
#   .[,.(nir_max5d = median(nir_max5d,na.rm=T), 
#        nir = median(nir,na.rm=T), 
#        nobs = sum(is.na(nir_max5d)==F)),
#     by=.(date,g)]
# ss[,`:=`(date = date(date))]
# ss <- ss[,`:=`(year=year(date),
#                month= month(date))]
# 
# 
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
#   geom_sf(data = amc, aes(fill = fzone)) +
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
# # ggsave(filename = "figures/AMC_zones.png",
# #        width = 11,
# #        height = 15,
# #        units = 'cm',
# #        dpi = 350)
# 
# 
# ## Map of Retrievals during October
# cdat[date %in% seq(ymd("2016-10-01"),ymd("2016-11-01"),by="1 day")] %>% 
#   ggplot(aes(x,y,fill=nir*scale_factor))+
#   geom_raster()+
#   coord_sf()+
#   scale_fill_viridis_c(option = 'B', 
#                        limits = c(0,0.2),
#                        oob = scales::squish) + 
#   labs(fill = "nir Reflectance", 
#        title = "MCD43A4 Adelaide Metro. Coastal Zones, October 2016", 
#        x = NULL,
#        y = NULL) +
#   facet_wrap(~date, 
#              nrow = 3) + 
#   theme_linedraw()+
#   theme(axis.ticks.x = element_blank(),
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         strip.text = element_text(size = 7),
#         panel.grid = element_blank(),
#         panel.background = element_rect(fill = 'grey30'), 
#         legend.position = 'bottom')
# ggsave(filename = "figures/MCD43A4_2016-10_nir.png",
#        width = 25,
#        height = 20,
#        units = 'cm',
#        dpi = 350)
# 
# 
# 
# merge(flow, ss, 
#       by=c("date","month","year")) %>% 
#   ggplot(aes(dv_max5d, nir_max5d, color = g))+
#   # geom_point()+
#   geom_smooth(se = F, 
#               method = 'lm') + 
#   facet_wrap(~year)
# 
# merge(flow[year==2016], ss, 
#       by=c("date","month","year")) %>% 
#   ggplot(aes(dv_max5d, nir, color = g))+
#   # geom_point()+
#   geom_smooth(se = F, 
#               method = 'lm')
# 
# 
# 
# (p1 <- flow[year==2016][month %between% c(8,12)] %>% 
#     ggplot(aes(date, dv_max5d))+
#     # geom_point()+ 
#     geom_line() + 
#     labs(x = NULL, 
#          y = "Rolling 5-day Max Daily Volume (ML)")+ 
#     scale_x_date(date_labels = "%Y-%m-%d") + 
#     theme_linedraw() + 
#     theme(panel.grid = element_blank()))
# 
# 
# 
# (p2 <- ss %>% 
#     ggplot(aes(date, nir_max5d*scale_factor,color=g))+
#     geom_line()+
#     geom_point(aes(size = nobs), 
#                alpha = 0.85, 
#                shape = 21) + 
#     scale_color_discrete_c4a_cat(palette = 'tol.muted') +
#     labs( y = "5-day Max nir Reflectance", 
#           x = NULL, 
#           size = "Number of pixel obs.",
#           color = "Zone and depth") + 
#     scale_x_date(date_labels = "%Y-%m-%d") +
#     theme_linedraw() + 
#     theme(panel.grid = element_blank(), 
#           legend.position = 'bottom'))
# 
# (p_out <- p1/p2)
# 
# ggsave(p_out,
#        filename = "figures/figure_eda_mcd43a4_",
#        width = 30,
#        height = 20,
#        units = 'cm',
#        dpi = 350)
# 
# 
# 
# 
# 
# 
