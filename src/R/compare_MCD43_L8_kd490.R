pacman::p_load(stars, lubridate, terra, tidyverse, 
               tidyterra,
               data.table,
               mgcv, mgcViz,sf, viridis, 
               cols4all, 
               patchwork)
getDTthreads()

# Torrens river discharge ------------------------
flow <- fread("data/daily_discharge_torrens_river_seaview_road_bridge.csv")
names(flow) <- c("date","time","dv")
flow[,`:=`(date = dmy(date))]
flow[,`:=`(year=year(date),month=month(date))]
flow <- flow[order(date)][,`:=`(dv_max5d = frollapply(dv,5,max,na.rm=T,align='right'))]
flow <- flow[order(date)][,`:=`(dv_max15d = frollapply(dv,15,max,na.rm=T,align='right'))]
flow <- flow[order(date)][,`:=`(dv_max30d = frollapply(dv,30,max,na.rm=T,align='right'))]

## Load and collate raster stacks of Landsat
vec_dates <- fread("data/gee_SAWater_turbidity/DATES_LC08_C02_T1_L2_kd490_l8_2013-03-18_2023-09-18.csv")
ic <- rast("data/gee_SAWater_turbidity/LC08_C02_T1_L2_kd490_l8_2013-03-18_2023-09-18.tif")

# ## extract depth and AMC zones ------------------
amc <- vect("data/AMC/AnalysisExtentMask_GDA94z54.shp")
amc <- terra::project(amc, crs(ic))
amc <- amc %>% select(zone)


time(ic) <- vec_dates$date %>% 
  as_date()
ic <- crop(ic,amc)
ic_l8 <- clamp(ic,lower=0,upper=20)



## Load and collate raster stacks of daily MODIS 
vec_dates <- fread("data/gee_SAWater_turbidity/DATES_MCD43A4_wang_kd490_daily_2000-02-24_2012-12-31_.csv")
vec_dates2 <- fread("data/gee_SAWater_turbidity/DATES_MCD43A4_wang_kd490_daily_2012-12-31_2023-07-11_.csv")

ic <- rast("data/gee_SAWater_turbidity/MCD43A4_wang_kd490_daily_2000-02-24_2012-12-31_.tif")
ic2 <- rast("data/gee_SAWater_turbidity/MCD43A4_wang_kd490_daily_2012-12-31_2023-07-11_.tif")

time(ic) <- vec_dates$date
time(ic2) <- vec_dates2$date
names(ic) <- vec_dates$date
names(ic2) <- vec_dates2$date

ic <- c(ic,ic2)
rm(ic2)

ic_mcd <- ic %>% clamp(.,lower=0,upper=20)


ss_mcd <- ic_mcd[[time(ic_mcd)
                  %between% 
            c(min(time(ic_l8)),max(time(ic_l8)))]]

ic_l8_c <- terra::project(ic_l8, ss_mcd, method='bilinear',threads=T)
ic_l8_c <- ic_l8_c[[time(ic_l8_c) %in% time(ss_mcd)]]


bath <- vect("data/Bathymetry/DetailedModel_DEPTH.shp")
names(bath) <- names(bath) %>% tolower()
r_bath <- rasterize(bath, ic_mcd,field="val_1")
names(r_bath) <- "depth"
r_amc <- rasterize(amc, ic_mcd,field="zone")
tmp0 <- c(r_bath,r_amc) %>% as.data.table(xy=T)

bath <- bath %>%
  mutate(depth_class = case_when(val_1 < 5 ~ "< 5",
                               val_1 >5 & val_1 <10 ~"5 - 10",
                               val_1 > 10 ~ "10"))
bath <- aggregate(bath,by='depth_class',dissolve=T)
bath <- bath %>% select(depth_class)
bath <- bath %>% filter(is.na(depth_class)==F)
crs(bath) <- crs(amc)



tmp1 <- ic_l8_c %>% 
  as.data.table(xy=T) %>%
  pivot_longer(-c("x","y")) %>% 
  mutate(time = ymd(str_extract(name,"\\d{8}"))) %>% 
  rename(kd_l8 = value) %>% 
  as.data.table()
tmp2 <- ss_mcd %>% 
  # stars::st_as_stars() %>% 
  as.data.table(xy=T) %>% 
  pivot_longer(-c("x","y")) %>% 
  rename(kd_mcd = value,
         time = name) %>% 
  as.data.table() %>% 
  mutate(time = as_date(time))

# microbenchmark::microbenchmark(as.Date("2023-01-01"))
# microbenchmark::microbenchmark(as_date("2023-01-01"))

dat <- merge(tmp1,tmp2,by=c('x','y','time'),all=T)
dat <- merge(dat, tmp0, by=c("x","y"),all=T)
dat <- dat %>% 
  mutate(depth_class = cut(depth,c(0,5,10,Inf)))

dat[,`:=`(zone = fct_relevel(zone, 
                  c("North Zone","Central Zone","South Zone")))] 


dat[,`:=`(id = .GRP), by=.(x,y)]
dat[,`:=`(year=year(time),
          month = month(time))]

pix_xy <- dat[,.(val = sum(is.na(kd_l8)==F)),by=.(x,y,id,zone,depth_class)][
  is.na(zone)==F][is.na(depth_class)==F][,
                  .(max_nobs = sum(val> 40)), by=.(zone,depth_class)]

dat <- merge(dat,pix_xy,by=c('zone','depth_class'),all=T)


# dat$kd_mcd %>% hist(100)
# dat$kd_l8 %>% hist(100)

tmp <- dat[is.na(zone)==F][depth>3][is.na(depth_class)==F] %>% 
  .[year %in% 2016:2017] %>% 
  .[,.(
  MCD = median(kd_mcd,na.rm=T),
  L8 = median(kd_l8,na.rm=T), 
  frac_nobs_l8 = sum(is.na(kd_l8)==F)/max_nobs),
                by=.(time,zone,depth_class)] %>% 
  pivot_longer(-c("time","frac_nobs_l8",'zone','depth_class')) %>% 
  setDT() %>% 
  .[is.na(value)==F]
  # .[frac_nobs_l8 > 0.33] %>% 


unique(tmp) %>% 
 ggplot(aes(time,value,color=name,group=name))+
  geom_line() + 
  geom_point(size=0.5) + 
  scale_color_brewer(type='qual',palette = 2) + 
  # scale_color_viridis_d(option='H') + 
  facet_grid(zone ~ depth_class,scales='fixed', 
             labeller = label_both) + 
  theme_linedraw() + 
  theme(axis.text.x = element_text(angle = 90))
ggsave(filename="figures/compare_mcd43_L8_kd490_by_depth.png",
       width=20,
       height=15,
       units='cm')

dat[is.na(kd_l8)==F] %>% 
  mutate(rmsd = sqrt((kd_mcd-kd_l8)**2)) %>% 
  ggplot(aes(depth,rmsd))+
  geom_point(data=. %>% sample_n(1000))+
  geom_smooth()

dat[is.na(kd_l8)==F] %>% 
  mutate(rmsd = sqrt((kd_mcd-kd_l8)**2)) %>% 
  ggplot(aes(kd_mcd,kd_l8,color=cut_number(depth,6)))+
  geom_smooth(method='lm',se=F)+
  geom_abline()+
  scale_color_viridis_d(option='H')


# dat[year %in% 2017][,.(
#   MCD = median(kd_mcd,na.rm=T),
#   L8 = median(kd_l8,na.rm=T), 
#   nobs_l8 = sum(is.na(kd_l8)==F)),
#   by=.(time,zone,depth_class)] %>% 
#   pivot_longer(-c("time","nobs_l8",'zone','depth_class')) %>% 
#   setDT() %>% 
#   .[is.na(zone)==F] %>% 
#   .[is.na(depth_class)==F] %>% 
#   # .[nobs_l8 > 2000] %>% 
#   # .[order(-name)] %>%
#   .[is.na(value)==F] %>% 
#   .[zone=="Central Zone" & depth_class == "(0,5]"] %>% 
#   .[name=="L8"]


dat[year==2016] %>% 
  ggplot(aes(kd_l8,kd_mcd))+
  geom_point() + 
  facet_wrap(~time)


dat[sample(.N,10000)] %>% 
  select(kd_mcd,kd_l8) %>% drop_na() %>% cor
  ggplot(aes(kd_mcd,kd_l8))+
  geom_point()+
  geom_smooth(method='lm') + 
  geom_abline()

ss1 <- dat[, .(l8 = median(kd_l8,na.rm=T),
        nobs_l8 = sum(is.na(kd_l8)==F),
        mcd = median(kd_mcd,na.rm=T)), 
    by=.(time)]
  
ss1 %>% 
  ggplot(aes(time, l8))+
  geom_point()+
  geom_line(aes(time,mcd))

ss1 %>% 
  ggplot(aes(time, mcd))+
  geom_line() + 
  geom_point(aes(time,l8,color=nobs_l8))+
  scale_color_viridis_c()


dat[, .(kd_8 = median(kd_l8,na.rm=T),
       kd_m = median(kd_mcd,na.rm=T)), 
    by=.(time)] %>% 
  .[m8<0.2] %>% 
  ggplot(aes(log10(m8),log10(t8)))+
  geom_point()+
  geom_smooth(method='lm')+
  geom_abline()



vec_ids <- sample(unique(dat$id), 30)

dat %>% 
  .[id %in% vec_ids] %>% 
  # [sample(.N,10000)][t_mcd<30][t_l8<30] %>% 
  # select(t_mcd, t_l8) %>% 
  # cor()
  ggplot(aes(t_mcd,t_l8, group=id))+
  geom_point() +
  geom_smooth(method='lm', se=F) + 
  geom_abline(col='#cf0000')




test <- dat[,`:=`(id = .GRP), by=.(x,y)]


s0 <- bam(t_l8 ~ s(month,bs='cc'), 
          data=dat, 
          select=T,
          discrete = T)
summary(s0); plot(s0,rug=T)
getViz(s0) %>% plot()
dat[is.na(t_l8)==F]$month %>% hist()
dat[is.na(t_mcd)==F]$month %>% hist()

s1 <- bam(t_mcd ~ s(month,bs='cc'), 
          data=dat, 
          select=T,
          discrete = T)
summary(s1); plot(s1)


s1 <- bam(t_l8 ~ 
            te(x,y,by=month), 
          data=dat %>% mutate(month= factor(month)), 
          select=T,
          discrete = T)
summary(s1)
getViz(s1) %>% plot()


m0 <- bam(log(kd_l8) ~ s(log(kd_mcd),bs='ts',k=5), 
          data=dat[sample(.N, 1e5)],
          select=T,
          discrete = T)
summary(m0); plot(m0)

m0 <- bam(log(kd_l8) ~
            kd_mcd + 
            log(kd_mcd), 
          data=dat[sample(.N, 1e5)],
          select=T,
          discrete = T)
summary(m0); plot(m0)
getViz(m0) %>% plot(allTerms=T) %>% print(pages=1)


m0 <- bam(log(kd_l8) ~
            s(depth) + 
            te(x,y)+ 
            s(month,bs ='cc') + 
            kd_mcd + 
            log(kd_mcd), 
          data=dat[sample(.N, 1e6)],
          select=T,
          discrete = T)
summary(m0)
getViz(m0) %>% plot(allTerms=T) %>% print(pages=1)


m0 <- bam(log10(kd_l8) ~
            te(x,y,month,bs=c("ts","ts","cc")) +
            # log10(kd_mcd) + 
            te(x,y,by=log10(kd_mcd))+
            # s(month,bs ='cc') +
            # kd_mcd +
            # log10(kd_mcd) +
            s(depth),
          data=dat[is.na(zone)==F][sample(.N, 1e6)][kd_l8>0.0001],
          select=T,
          discrete = T)
summary(m0)
getViz(m0) %>% plot(allTerms=T) %>% print(pages=1)

m1 <- gam(log10(kd_l8) ~
            te(x,y,month,bs=c("ts","ts","cc")) +
            # log10(kd_mcd) + 
            log10(kd_mcd)+ 
            s(x,y,by=log10(kd_mcd))+
            # s(month,bs ='cc') +
            # kd_mcd +
            # log10(kd_mcd) +
            s(depth),
          data=dat[is.na(zone)==F][sample(.N, 5e5)][kd_l8>0.0001][depth>1],
          select=T,
          # discrete = T
          method='REML'
          )
summary(m1)
getViz(m1) %>% plot(allTerms=T) %>% print(pages=1)



# dat$kd_l8 %>% log() %>% hist(100)
# 
# dat[sample(.N,10000)][t_mcd<30][t_l8<30] %>% 
#   select(t_mcd, t_l8) %>% 
#   cor()
#   ggplot(aes(t_mcd,t_l8))+
#   geom_point() +
#   geom_smooth(method='lm') + 
#   geom_abline(col='#cf0000')
# 
# 
# # tmp <- c(ic_l8_c,ss_mcd) %>% 
# #   stars::st_as_stars() %>% 
# #   as.data.table(xy=T)
# # 
# # tmp
# 
# dim(ss_mcd)
# dim(ic_l8_c)
# 
# 
# 
# 
# 
# 
# # 
# # bath <- vect("data/Bathymetry/DetailedModel_DEPTH.shp")
# # names(bath) <- names(bath) %>% tolower()
# # bath <- bath %>% 
# #   mutate(depth_class = case_when(val_1 < 5 ~ "< 5",
# #                                val_1 >5 & val_1 <10 ~"5 - 10",
# #                                val_1 > 10 ~ "10"))
# # bath <- aggregate(bath,by='depth_class',dissolve=T)
# # bath <- bath %>% select(depth_class)
# # bath <- bath %>% filter(is.na(depth_class)==F)
# # crs(bath) <- crs(amc)
# # 
# # bath
# # 
# # tmp
# # aggregate(merge(amc,bath),by=c("zone","depth_class"))
# # 
# # vect(union(amc,bath))
# # vect(intersect(amc,bath))
# # vect(intersect(bath,amc))
# # 
# # union(amc,bath)[[1]]$depth_class
# # union(amc,bath)[[2]]$depth_class
# # union(amc,bath)[[3]]$depth_class
# # union(amc,bath)[[4]]$depth_class
# # union(amc,bath)[[5]]$depth_class
# # union(amc,bath)[[6]]$depth_class
# # union(amc,bath)[[7]]$depth_class
# 
# # tmp <- union(amc,bath) %>% 
# #   vect %>% 
# #   aggregate(., by=c("zone","depth_class"))
# # tmp %>% filter(is.na(zone)==F) %>% 
# #   filter(is.na(depth_class)==F)
# # 
# # aggregate(tmp) %>% plot
# # 
# # plot(tmp)
# # names(tmp)
# # 
# # 
# # grid <- rast("data/gee_SAWater_turbidity/LC08_C02_T1_L2_turb_l8_2013-03-18_2023-09-18.tif")[[1]]
# # bath <- rasterize(bath, grid, field="depth_class")
# # 
# # amc <- rasterize(amc,grid,field='zone')
# # 
# # tmp <- c(bath,amc)
# # terra::crosstab(tmp)
# # 
# # 
# # 
# # as.data.table(xy=T) %>% 
# #   .[is.na(zone)==F]
# # tmp %>% 
# #   as.polygons()
# #   plot()
# 
# 
# 
# 
# 
# ic[[year(time(ic))==2016]] %>% 
#   plot()
# 
# # tmp <- ic[[time(ic)==ymd("2016-12-08")]]
# # 
# # tmp2 <- quantize_raster(tmp %>% log10(), 
# #                 n_levels = 16, method="equal prob") 
# # 
# # tmp2 %>% 
# #   plot(col=grey.colors(16))
# # freq(tmp2)[,c("value", "count")]
# # 
# # textures1 <- glcm_textures(tmp2, 
# #                           w = c(3,3), 
# #                           metrics = "glcm_entropy",
# #                           n_levels = 16, 
# #                           quantization = "none", 
# #                           shift = c(1,0)) 
# # plot(textures1,col=viridis::viridis(10))
# # library(fastGLCM)
# 
# 
# v <- 2013:2023 %>% 
#   lapply(., FUN= function(x){
#     out <- ic[[x == year(time(ic))]] %>% stdev(na.rm=T)
#     time(out) <- x
#     return(out)
#   })
# 
# 
# 
# pal <- inferno(150)
# plot(ic[[time(ic)==ymd("2016-10-05")]] %>% clamp(., lower=0,upper=50),
#     colNA='grey',
#     maxcell=1e6,
#     smooth=T,
#    col=pal,
#    type='continuous',
#    range=c(0,30)
#    )
# 
# tmp <- v[[4]] %>% focal(., w=5, fun='median')
# plot(v[[4]],col=pal,range=c(0,10))
# plot(tmp,col=pal,range=c(0,15),colNA='grey30')
# 
# 
# plot(rast(v),col=pal,range=c(0,10))
# 
# 
# 
# ggplot()+
#   geom_spatraster(data=ic[[time(ic)==ymd("2016-10-05")]])+
# scale_fill_viridis_c(option='C')
# 
# 
# 
# vec_dates$date
# 
# 
# 
# 
# 
# ic_50 <- median(ic,na.rm=T)
# ic_sd <- stdev(ic,na.rm=T)
# 
# plot(ic_sd,col=pal,colNA='grey30')
# 
# flow[year==2016][month%in%c(9,10)] %>% 
#   ggplot(aes(date,dv))+
#   geom_line()
# 
# 
# 
# # ic2 <- rast("data/gee_SAWater_turbidity/LC08_C02_T1_L2_turb_l8_2013-03-18_2023-09-18-0000002048-0000000000.tif")
# # 
# # ic <- terra::merge(ic1,ic2)
# # terra::mosaic(ic,ic2)
# 
# 
# 
# 
# 
# amc2 <- union(bath,amc) %>% vect
# amc2 <- amc2[is.na(amc2$'zone')==F] 
# 
# plot(amc2)
# 
# aggregate(amc2,by)
# 
# aggregate(bath,amc)
# 
# zonal(bath,amc,fun='unique')
# 
# 
# st_combine(bath)
# 
# 
# bath <- terra::project(bath, crs(amc))
# 
# 
# 
# ic_u <- mean(ic,na.rm=T)
# plot(ic_u,col=pal,colNA='grey30')
# 
# ic_u %>% summary
# 
# 
# ic_max <- max(ic,na.rm=T)
# 
# ## Define coords
# xy <- ic_u %>% 
#   as.data.table(xy=T) %>% 
#   select(x,y)
# 
# 
# 
# 
# bath <- rasterize(bath, ic, field="Val_1")
# names(bath) <- "depth"
# 
# # amc <- rasterize(amc, ic, field="zone")
# 
# 
# 
# tmp0 <- merge(bath, amc) %>% 
#   mutate(depth_class = case_when(Val_1 < 5 ~ "< 5",
#                                  Val_1 >5 & Val_1 <10 ~"5 - 10",
#                                  Val_1 > 10 ~ "10"))
# aggregate(tmp0 %>% select("depth_class","zone"), 
#           by=c("depth_class","zone"),dissolve=TRUE) %>% 
#   plot()
# tmp0
# 
# 
# tmp <- merge(bath, amc) %>% 
#   st_as_sf() %>% 
#   mutate(depth_class = case_when(Val_1 < 5 ~ "< 5",
#                                  Val_1 >5 & Val_1 <10 ~"5 - 10",
#                                  Val_1 > 10 ~ "10"))
# tmp %>% 
#   select(depth_class,zone) %>% 
#   plot()
# 
#     
# aggregate(.,
#             # fun = 'mean',
#             by=c('depth_class','zone'))
# 
# bath %>% 
#   # st_as_sf() %>% 
#   mutate(depth_class = case_when(Val_1 < 5 ~ "< 5",
#                                  Val_1 >5 & Val_1 <10 ~"5 - 10",
#                                  Val_1 > 10 ~ "10")) %>% 
#   select(depth_class) %>% 
#   aggregate(., by='depth_class') %>% 
#   plot()
# 
# b5 <- subset(bath,bath$"Val_1" < 5)
# b5_10 <- subset(bath,bath$"Val_1" < 10 & bath$"Val_1">5)
# b10 <- subset(bath,bath$"Val_1" > 10)
# 
# merge(b5 %>% mutate(dclass = "< 5"),
# b5_10 %>% mutate(dclass = "5 - 10")) %>% plot
# 
# c(b5 %>% mutate(dclass = "< 5"),
# b5_10 %>% mutate(dclass = "5 - 10"),
# b10 %>% mutate(dclass = "> 10"))
# 
# vv_bath <- terra::extract(bath, xy[,.(x,y)], xy=T) %>% 
#   as.data.table()
# vv_amc <- terra::extract(amc, xy[,.(x,y)], xy=T) %>% 
#   as.data.table()
# 
# # 
# # tmp <- terra::extract(ic, amc, fun=max,na.rm=T)
# # tmp %>% 
# #   pivot_longer(-ID) %>% 
# #   mutate(time = ymd(str_extract(name, "[0-9]{8}"))) %>% 
# #   pull(value)
# # 
# # zonal(amc, ic)
# # zonal(ic,amc,"mean",na.rm=T)
# # 
# 
# 
# tmp <- ic %>% 
#   as.data.table(xy=T)
# tmp2 <- tmp %>% 
#   pivot_longer(-c("x","y"))
# setDT(tmp2)
# 
# tmp3 <- tmp2 %>% as.data.table() %>% .[,.(val50 = median(value,na.rm=T)),by=.(x,y)]
# tmp4 <- rast(tmp3,type='xyz')
# 
# ggplot()+
#   geom_spatraster(data=tmp4)+
#   geom_sf(data=st_as_sf(amc), 
#           fill='transparent')+
#   scale_fill_viridis_c(option='B',
#                        limits=c(0,30))+
#   coord_sf()
# 
# tmp_nobs <- tmp2[,.(nobs = sum(is.na(value)==F)),by=.(x,y)]
# ggplot()+
#   geom_spatraster(data=rast(tmp_nobs,type='xyz'))+
#   geom_sf(data=st_as_sf(amc), 
#           fill='transparent')+
#   scale_fill_viridis_c(option='B',
#                        limits=c(0,30))+
#   coord_sf()
# 
# 
# tmp_nobs[sample(.N,50000)] %>% 
#   ggplot(aes(x,y,color=nobs))+
#   geom_point(size=0.3)+
#   scale_color_viridis_c(option='B',
#                         limits=c(60,120)) + 
#   coord_sf()
# 
# 
# 
# ic <- clamp(ic, lower = 0, upper = 100)
# plot(ic[[1]])
# d5 <- exact_extract(ic, st_as_sf(amc), fun = c('mean','count','median'), 
#                     append_cols = "zone")
# 
# 
# d5 %>% 
#   pivot_longer(-zone) %>% 
#   as.data.table() %>% 
#   mutate(time = ymd(str_extract(name, "\\d{8}"))) %>% 
#   rowwise() %>% 
#   mutate(stat = str_split(name,"\\.",simplify = T)[1]) %>% 
#   as.data.table() %>% 
#   filter(stat == 'median') %>% 
#   ggplot(aes(time,value,color=zone))+
#   geom_line()
