pacman::p_load(stars, sf, terra, data.table, tidyverse, lubridate,
               arrow, mgcv, gratia, mgcViz, patchwork, 
               phenofit)

# NOTES: 
"
* Port River dredging in early 2006
* Penrice peak discharge was in 2011 and stayed high until closing in 2013
* Two sludge outfalls closed in 1993
* Glenelg sludge discharge finished in 2011
* drought 2001-2009, followed by increased suspended solids in 2010
* 
* Questions: 
* What occurred in late 2019/2020? Could the turbidity peak be 
  connected to the fires?
"


# Define Functions   ==============================================
## Dogliotti turbidity
fn_dogliotti <- function(p_w){
  # Dogliotti et al 2015, Remote Sensing of Environment
  # Notes: See ACOLITE for most recent parameter calibrations
  A_t <- 378.46
  C <- 0.19905
  out <- A_t * (p_w) / 
    (1 - ((p_w)/ C))
  return(out)
}

fn_nechad <- function(p_w){
  # Nechad et al., 2010 
  # Notes: See parameter table in Nechand 2010 for sensor specific Cp
  ## Table 2. Cp: MODIS "HIRES" 645: 16.41
  ## Table 3. Bp: 645nm: 2.32
  ## Table 3. Ap: 645nm: 253.51
  Cp <- 16.41
  Bp <- 2.32 
  Ap <- 253.51
  S <- (Ap*p_w)/(1 - (p_w/Cp)) + Bp
  return(S)
}


# Import Data  ==============================================
dat <- arrow::read_parquet("outputs/MCD43A4_AMC_red.parquet") %>% 
  mutate(year = year(time),
         month = month(time), 
         doy = yday(time)) %>% 
  mutate(fzone = factor(dz)) %>% 
  mutate(dog_mean = fn_dogliotti(mean),
         nechad_mean = fn_nechad(mean))



# Smooth time series ==========================================
vec_zones <- dat$fzone %>% levels()
d_time <- data.table(time = seq(dat[mean>0]$time %>% min,
                      to = dat[mean>0]$time %>% max,
                      by='1 day'))

fn_smooth <- function(x){
  tmp_dat <- dat[fzone==x]
  tmp_dat <- merge(tmp_dat, d_time,
                   by= 'time',
                   # allow.cartesian=T, 
                   all = T) %>% 
    mutate(year=year(time),
           month=month(time),
           doy = yday(time))

  tmp_dat <- tmp_dat[order(time)]
  tmp_dat <- tmp_dat %>% 
    mutate(w = count) %>% 
    mutate(w = ifelse(is.na(w)==T, 0, w))
  

  tmp_dat[dog_mean==0]$dog_mean <- NA_real_
  y <- tmp_dat$dog_mean
  y <- nafill(y,type='locf')
  y <- nafill(y,type='nocb')
  y <- nafill(y,type='locf')
  ylu <- range(y)
  
  sy <- phenofit::smooth_wWHIT(y,
                         w = tmp_dat$w, 
                         ylu = ylu,
                         nptperyear = 365, 
                         iters = 10, 
                         lambda = NULL)
  tmp_dat <- tmp_dat %>% mutate(sdog_mean = sy$zs$ziter10) %>% 
    mutate(fzone = x)
  return(tmp_dat)
}

# sy$zs$ziter10
# 
# # sy$zs$ziter1 %>% plot(type='l')
# # y %>% points
# # sy$zs$ziter10 %>% lines(col='red')
# x <- vec_zones[1]
# 
# fn_smooth(vec_zones[1]) %>% 
#   pull(smean)
d2 <- vec_zones %>% lapply(., fn_smooth)
d2 <- rbindlist(d2) %>% 
  mutate(fzone = factor(fzone))

# fn_smooth('CentralZone_lt10')[year==2000 & month==4]
# d2[fzone=='CentralZone_lt10'][year==2000 & month==4]


# dat[fzone=='CentralZone_lt10'] %>% dim
# fn_smooth('CentralZone_lt10') %>% dim
# fn_smooth('CentralZone_lt10')[time >= ymd("2000-03-27")][time<=ymd("2000-04-30")]$sdog_mean
# 
# dat[fzone=='CentralZone_lt10'][time >= ymd("2000-03-27")][time<=ymd("2000-04-30")]$dog_mean
# d2[fzone=='CentralZone_lt10'][time >= ymd("2000-03-27")][time<=ymd("2000-04-30")]$dog_mean
# d2[fzone=='CentralZone_gt10'][time >= ymd("2001-03-27")][time<=ymd("2001-05-01")]$sdog_mean
# d2$fzone %>% table


# calculate anomaly ===========================================
m0 <- bam(log(sdog_mean) ~ 
            s(doy,fzone,bs=c("fs","cc"),k=5),
            # s(doy,by=fzone,bs='cc'), 
          data=d2, 
          weights = d2$count,
          select = F, 
          discrete = T)
summary(m0)
# plot(m0)
# draw(m0)

gratia::smooth_estimates(m0) %>% 
  ggplot(aes(doy, .estimate,color=fzone))+
  geom_line()

ref <- expand.grid(doy = 1:365,
                   fzone = levels(dat$fzone)) %>% setDT()
ref <- ref %>% mutate(sdog_u = predict(m0, newdata=., type='response') %>% 
                        exp())

d3 <- merge(d2,ref,by=c("fzone","doy"))
d3[,`:=`(sdog_mean_anom = sdog_mean - sdog_u)]



## Import daily stream discharges ------------------------
proc_flow <- function(fnp){
  tmp <- fread(fnp)
  names(tmp) <- c("date","time","dv")
  tmp[,`:=`(time = dmy(date))]
  tmp[,`:=`(year=year(time),month=month(time),doy=yday(time))]
  tmp <- tmp[order(time)][,`:=`(dv_max5d = frollapply(dv,5,max,na.rm=T,align='right'))]
  tmp <- tmp[order(time)][,`:=`(dv_max10d = frollapply(dv,10,max,na.rm=T,align='right'))]
  tmp <- tmp[order(time)][,`:=`(dv_max15d = frollapply(dv,15,max,na.rm=T,align='right'))]
  tmp <- tmp[order(time)][,`:=`(dv_max20d = frollapply(dv,20,max,na.rm=T,align='right'))]
  tmp <- tmp[order(time)][,`:=`(dv_max30d = frollapply(dv,30,max,na.rm=T,align='right'))]
  tmp <- tmp[order(time)][,`:=`(dv_max40d = frollapply(dv,40,max,na.rm=T,align='right'))]
  tmp <- tmp[order(time)][,`:=`(dv_max45d = frollapply(dv,45,max,na.rm=T,align='right'))]
  tmp <- tmp[order(time)][,`:=`(dv_max60d = frollapply(dv,60,max,na.rm=T,align='right'))]
  tmp <- tmp[order(time)][,`:=`(dv_max90d = frollapply(dv,90,max,na.rm=T,align='right'))]
  
  tmp <- tmp[order(time)][,`:=`(dv_sum5d = frollapply(dv,5,sum,na.rm=T,align='right'))]
  tmp <- tmp[order(time)][,`:=`(dv_sum10d = frollapply(dv,10,sum,na.rm=T,align='right'))]
  tmp <- tmp[order(time)][,`:=`(dv_sum15d = frollapply(dv,15,sum,na.rm=T,align='right'))]
  tmp <- tmp[order(time)][,`:=`(dv_sum20d = frollapply(dv,20,sum,na.rm=T,align='right'))]
  tmp <- tmp[order(time)][,`:=`(dv_sum30d = frollapply(dv,30,sum,na.rm=T,align='right'))]
  tmp <- tmp[order(time)][,`:=`(dv_sum40d = frollapply(dv,40,sum,na.rm=T,align='right'))]
  tmp <- tmp[order(time)][,`:=`(dv_sum45d = frollapply(dv,45,sum,na.rm=T,align='right'))]
  tmp <- tmp[order(time)][,`:=`(dv_sum60d = frollapply(dv,60,sum,na.rm=T,align='right'))]
  tmp <- tmp[order(time)][,`:=`(dv_sum90d = frollapply(dv,90,sum,na.rm=T,align='right'))]
  tmp <- tmp[order(time)][,`:=`(diff1 = dv - lag(dv,n = 1))]
  tmp <- tmp %>% select(-date)
  tmp[,`:=`(q1 = quantile(dv,0.1,na.rm=T))]
  return(tmp)
}



list.files("data","daily",full.names = T)
dflow <- rbindlist(list(
  proc_flow("data/daily_discharge_torrens_river_seaview_road_bridge.csv") %>% 
    mutate(site = "torrens_seaview"),
  proc_flow("data/daily_discharge_torrens_river_holbrooks_road.CSV") %>% 
    mutate(site = "torrens_holbrook"),
  proc_flow("data/daily_discharge_helps_drain_summer_bolivar.csv") %>% 
    mutate(site = "helps_drain"),
  proc_flow("data/daily_discharge_gawler_river_virginia_park.csv") %>% 
    mutate(site = "gawler_virginia"), 
  proc_flow("data/daily_discharge_christie_creek.csv") %>% 
    mutate(site = "christie_creek")
  
))



# d3[fzone=='CentralZone_lt10'][time >= ymd("2000-03-27")]$sdog_mean %>% 
#   is.na() %>% 
#   table()

tmp_dat <- merge(d3[fzone=='CentralZone_lt10'], 
              dflow[site=="torrens_holbrook"], 
              by=c("time","year","month","doy"), 
              all=T, 
              allow.cartesian = T)
# tmp_dat <- tmp_dat[time>=ymd("2000-03-27")]
# tmp_dat$sdog_mean %>% is.na() %>% table()
# tmp_dat[is.na(sdog_mean)==T]
# d3[fzone=='CentralZone_lt10'][time%in%tmp_dat[is.na(sdog_mean)==T]$time]
# d3[fzone=='CentralZone_lt10'][time==ymd("2000-04-03")]
# d3[fzone=='CentralZone_lt10'][year==2000 & month==4]


# Set up timeseries =================================
tmp_dat <- merge(d3[fzone=='CentralZone_lt10'], 
                 dflow[site=="torrens_holbrook"], 
                 by=c("time","year","month","doy"), 
                 all=T, 
                 allow.cartesian = T) %>% 
  filter(is.na(sdog_mean)==F) %>% 
  filter(time <= ymd("2023-08-09"))



ts1 <- ts(data = tmp_dat$dv, 
          frequency = 365, 
          start = c(2000,87))
ts2 <- ts(data = tmp_dat$sdog_mean, 
          frequency = 365, 
          start = c(2000,87))
ts3 <- ts(data = tmp_dat$sdog_mean_anom, 
          frequency = 365, 
          start = c(2000,87))

png(filename = 'figures/decompose_torrens_holbrook_dv.png',
    width = 15,
    height = 15,
    units='cm',
    pointsize = 10,
    res = 350)
decompose(ts1,type='additive') %>% 
  plot()
dev.off()

png(filename = 'figures/decompose_mcd43a4_smoothed_dogliotti_CentralZone_lt10.png',
    width = 15,
    height = 15,
    units='cm',
    pointsize = 10,
    res = 350)
decompose(ts2,type='additive') %>% 
  plot()
dev.off()


trend_dv <- decompose(ts1,type='additive')$trend
trend_sdog <- decompose(ts2,type='additive')$trend
trend_sdog_anom <- decompose(ts3,type='additive')$trend

plot(trend_dv,ylim=c(0,250))
lines(trend_sdog_anom*50,col='red')


vv1 <- data.table(ddate = time(ts1) %>% as.numeric(), 
                  trend_dv = trend_dv,
                  trend_sdog = trend_sdog,
                  trend_sdog_anom = trend_sdog_anom)
vv1 <- vv1[order(ddate)]

corr <- function(y) cor(y[, "trend_dv"], y[, "trend_sdog_anom"])
vv1[,rcor:=zoo::rollapplyr(.SD, 365, corr, by.column = FALSE, fill = NA)]
vv1$rcor %>% plot


decimal_date(ymd("2020-8-01"))

vv1[is.na(trend_dv)==F][ddate >= 2002][,.(ddate,trend_dv,trend_sdog_anom,rcor)] %>% 
  set_names("ddate","Torrens Discharge","Dogliotti FNU Anomaly","365 day Moving Correlation") %>% 
  pivot_longer(-ddate) %>% 
  ggplot(aes(ddate,value))+
  geom_line() + 
  geom_vline(aes(xintercept = decimal_date(ymd("2019-06-01"))), 
             col='red',lty=3) + 
  geom_vline(aes(xintercept = decimal_date(ymd("2019-08-01"))), 
             col='red',lty=3) + 
  geom_hline(aes(yintercept = 0),
             lty=2)+
  labs(title = "Central Zone, < 10 m") + 
  facet_wrap(~name, ncol = 1, 
             scales = 'free') +
  theme_linedraw()+
  theme(panel.grid.minor = element_blank())
ggsave(filename = "figures/fig_timeseries_dogAnom_TorrensDischarge_centralZone_lt10.png",
       width = 20,
       height = 15,
       units='cm')


# Autoregressive models =================================
ts_dv <- ts(data = tmp_dat$dv, 
          frequency = 365, 
          deltat = 1/365,
          start = c(2000,87))
ts_sdog <- ts(data = tmp_dat$sdog_mean, 
          frequency = 365, 
          start = c(2000,87))
ts_sdog_anom <- ts(data = tmp_dat$sdog_mean_anom, 
          frequency = 365, 
          start = c(2000,87))


ts_dog <- ts(data = tmp_dat$dog_mean, 
             frequency = 365, 
             start = c(2000,87))
ar(ts_dog)

w1 <- arima(ts_dog)
w1

ts_dog %>% plot
acf(ts_dog)
acf(ts_sdog_anom)


d_tmp2
mm1 <- gamm(
  trend_sdog_anom ~ s(trend_dv,k=5),
  correlation = corAR1(),
 data = d_tmp2[order(ddate)]
)
mm1$lme
mm1$lme
plot(mm1$gam)
acf(residuals(mm1$gam),main="raw residual ACF")
acf(residuals(mm1$lme,type="normalized"),main="standardized residual ACF")

plot(gam(trend_sdog_anom ~ s(trend_dv,k=5),
         data=d_tmp2))




tmp_dat[,`:=`(ddate = decimal_date(time))]
d_tmp3 <- tmp_dat[time >= ymd("2002-06-01")][time<=ymd("2003-06-01")]

head(d_tmp3)
mm1 <- gamm(
  sdog_mean_anom ~ dv,
  correlation = corAR1(),
  data = d_tmp3[order(time)]
)
summary(mm1$gam); 


glmmTMB(sdog_mean_anom ~ dv + ar1(ddate + 0|site), 
        data=tmp_dat[1:100,])


d_tmp3$dog_mean
gls(dog_mean ~ dv, 
    correlation = corARMA(c(0.5,0.5), p=1,q=1), 
    data=d_tmp3)

names(d_tmp3)
gamm(dog_mean ~ s(dv,k=5), 
     # correlation = corAR1(form = ~ddate),
      correlation = corARMA(form = ~ddate|year),
     data=d_tmp3)


decompose(ts3,type='additive') %>% 
  plot()


tmp_dat$doy[1]
decompose(ts1,type='multiplicative') %>% plot()
stl(ts1, s.window = "periodic") %>% plot
# StructTS(ts1)
# ?plot.decomposed.ts

decompose(ts2) %>% plot


# Granger Causality ======================
tsDat <- ts.union(ts1, ts2) 
tsDat
tsVAR <- vars::VAR(tsDat, p = 1, lag.max=45)
vars::causality(tsVAR, cause = "ts1")$Granger




# ts_4 <- ts(rnorm(1000))
# ts_5 <- ts(rnorm(1000))
# vars::VAR(ts.union(ts_4,
#                    ts_5), 
#           p=1) %>% 
#   vars::causality(cause = 'ts_4')
# 
# 
# lmtest::grangertest(ts1, ts2, order=1)
# lmtest::grangertest(ts1, ts2, order=30)
# 
# ccf(ts1,ts2,lag.max=90,type='correlation')
# ccf(ts1,ts2,lag.max=90,type='covariance')
# ccf(ts1,ts2,lag.max=90,type="covariance", 
#     plot=F,
#     demean=T)
# 
# 
# 
# 
# tmp_dat[is.na(sdog_mean)==F]$time %>% summary
# tmp_dat[time >= ymd("2000-03-27")][time<=ymd("2000-04-30")]$sdog_mean
# 
# d2[time >= ymd("2000-03-27")][time<=ymd("2000-04-30")]$sdog_mean
# 
# 
# 
# tmp_dat %>% 
#   ggplot(aes(doy, sdog_mean_anom,color=year,group=year))+
#   geom_line() +
#   geom_smooth(aes(color=NULL,group=NULL), 
#               method='lm') + 
#   scale_color_viridis_c()
# 
# 


# Corrplot w/different streamflow lags ================
sfac <- 2
png(filename = "figures/corrplot_mcd43a4_dog_flow-sum.png",
    width = 10*sfac,
    height = 10*sfac, 
    units = 'cm',
    pointsize = 10,
    res = 300)
tmp_dat %>% 
  rename(fnu_mean = sdog_mean, 
         fnu_mean_anom = sdog_mean_anom) %>% 
  select(fnu_mean,fnu_mean_anom, starts_with("dv_sum")) %>% 
  drop_na() %>% 
  cor() %>% 
  corrplot::corrplot(method='number', 
                     type = 'upper', 
                     title = "MCD43A4 Dogliotti Turbidity, Central Zone < 10 m", 
                     mar = c(0,0,1,0),
                     diag = F)
dev.off()

sfac <- 2
png(filename = "figures/corrplot_mcd43a4_dog_flow-max.png",
    width = 10*sfac,
    height = 10*sfac, 
    units = 'cm',
    pointsize = 10,
    res = 300)
tmp_dat %>% 
  rename(fnu_mean = sdog_mean, 
         fnu_mean_anom = sdog_mean_anom) %>% 
  select(fnu_mean,fnu_mean_anom, starts_with("dv_max")) %>% 
  drop_na() %>% 
  cor() %>% 
  corrplot::corrplot(method='number', 
                     type = 'upper', 
                     title = "MCD43A4 Dogliotti Turbidity, Central Zone < 10 m", 
                     mar = c(0,0,1,0),
                     diag = F)
dev.off()



# Holt Winters Exponential Smoothing
## s[0] = x[0]
## s[t] = alpha*x[t] + (1-alpha)*s[t-1], t > 0
## alpha is the exponential smoothing factor, 0 < alpha < 1

m2 <- HoltWinters(ts2, gamma=F, beta=F)
(m2)
plot(m2)


?forecast::ets
plot(ts2)
m3 <- forecast::ets(ts1, 
    model = "ANN")
summary(m3)
plot(m3)

m4 <- forecast::ets(ts2, 
                    model = "ZZZ")
m4

m5 <- ar(ts2)
m5
predict(m5,n.ahead=365)$pred %>% plot


?tsDyn::nlar

tsDyn::linear(lynx,m=3)
tsDyn::linear(lynx,m=2, type='ADF') %>% plot
mod.nnet <- nnetTs(log(lynx), m=2, size=3)
mod.nnet

# tmp_dat %>% 
#   bam(sdog_mean_anom ~ 
#         s(doy) + 
#         # s(dv_sum5d) + 
#         # s(dv_max5d) + 
#         # s(dv_sum45d, k=5) +
#         s(diff1) + 
#         s(I(dv_max45d/dv_sum45d), k=5),
#       data=.,
#       select = T,
#       discrete = T) %>% 
#   # summary()
#   plot(scheme=2, pages=1)
# 
# tmp_dat[sdog_mean > 0]
# 
# 
# 
# 
# 
# 
# d3$count %>% hist(100)
# d3[count > 10] %>% 
#   ggplot(aes(time, dog_mean))+
#   geom_point() + 
#   geom_smooth(formula = y~s(x,bs='ad',k=50),
#               aes(weight = count))
# 
# d3[dv_sum45d>0][year >= 2000] %>% 
#   ggplot(aes(time, dv_sum45d))+
#   geom_point() + 
#   geom_smooth(formula = y~s(x,bs='ad',k=50), 
#               method = 'gam',
#               method.args = list(family = Gamma(link='log')))
# 
# 
# 
# d3[dv_sum45d==0]
# b1 <- gam(dog_mean_anom ~ 
#             s(val,k=5),
#           data=d3 %>% 
#             filter(dv_sum45d > 0) %>% 
#             mutate(val = (dv_sum45d)),
#           method='REML')
# summary(b1); plot(b1,pages=1)
# 
# 
# d3[year%in%c(2019,2020)] %>% 
#   plot(dog_mean_anom~time, data=.,type='l')
# d3[year%in%c(2019,2020)] %>% 
#   lines(dv_sum45d/300~time, data=.,type='l',col='red')
# 
# 
# 
# d3[year==2020] %>% 
#   ggplot(aes(dv_sum40d, mean,color=doy))+
#   geom_point()+
#   geom_path() + 
#   scale_color_viridis_c(option='H')
# 
# 
# d3 %>% 
#   select(mean,starts_with("dv_max")) %>% 
#   drop_na() %>% 
#   cor() %>% 
#   corrplot::corrplot(method='number')
# 
# 
# dflow %>% 
#   ggplot(aes(time, dv_max5d, color=site))+
#   geom_line()+
#   facet_wrap(~site)
# 
# 
# 
# 
# 
# tmp <- dflow[site=='torrens_holbrook'][order(time)][,`:=`(
#   dv_l1 = shift(dv,n = 1)
# )][,`:=`(flow_presence = case_when(
#   dv > 0 ~ T,
#   .default = F
# ))]
# tmp[,`:=`(flow_start = case_when(
#   (dv==dv_max5d)&(dv>0.5)&(flow_presence==T) ~ T,
#   .default = F
# ))]
# tmp[flow_start==T]
# tmp[flow_start==T][year==2016] %>% 
#   ggplot(aes(time, dv))+
#   geom_line(data=tmp[year==2016],
#             aes(time,dv)) + 
#   geom_vline(aes(xintercept = time),
#              col='red')
# 
# dat2 <- merge(dat, 
#                dflow[site=="torrens_holbrook"], 
#                by=c("time","year","month","doy"), 
#               all=T, 
#               allow.cartesian = T)
# dat2 <- full_join(dat, 
#               dflow[site=="torrens_holbrook"], 
#               by=c("time","year","month","doy"))
# 
# 
# dflow[year==2020 &site=="torrens_holbrook"]$dv_sum15d
# dat2[year==2020 &site=="torrens_holbrook"&fzone=="CentralZone_lt10"]$dv_sum15d
# 
# (dat2[fzone=="CentralZone_lt10"][year==2020] %>% 
#   ggplot(aes(time, dog_mean_anom))+
#   geom_point()) /
#   dflow[year==2020 &site=="torrens_holbrook"] %>% 
#   ggplot(aes(time, dv_sum30d))+
#   geom_line()
# 
# c1 <- bam(dog_mean_anom ~ s(dv_sum90d,k=10), 
#           data=dat2[fzone=="CentralZone_lt10"],
#           discrete= T,
#           select=T)
# summary(c1)
# plot(c1)
# 
# 
# tmp
# tmp[flow_presence==T]
# tmp[year==2020] %>% 
#   ggplot(aes(time, dv))+
#   geom_line() + 
#   geom_vline(data=. %>% filter(flow_start==T),
#              aes(xintercept = time))
# tmp[flow_start==T] %>% 
#   ggplot(aes(time, dv))+
#   geom_line() + 
#   geom_vline(data=. %>% filter(flow_start==T),
#              aes(xintercept = time),
#              col='red')
# 
# 
# 
# # dcast(dflow, site+time+year+month+doy ~ dv, fun=mean)
# # dcast(dflow, site+time ~ dv, fun=mean)
# # pivot_wider(dflow, 
# #             names_from = 'site',
# #             id_cols = c("time"),
# #             values_from = 'dv') %>% 
# #   ggplot(aes(gawler_virginia, torrens_holbrook))+
# #   geom_point()
# 
# m0 <- bam(mean ~ 
#             # s(month,fzone,bs=c("fs","cc")) + 
#             s(doy,by=fzone,bs='cc'), 
#           data=dat, 
#           weights = dat$count,
#           select = F, 
#           discrete = T)
# summary(m0)
# plot(m0)
# draw(m0)
# 
# gratia::smooth_estimates(m0) %>% 
#   ggplot(aes(doy, .estimate,color=fzone))+
#   geom_line()
# 
# ref <- expand.grid(doy = 1:365,
#            fzone = levels(dat$fzone)) %>% setDT()
# ref <- ref %>% mutate(dog_u = predict(m0, newdata=., type='response'))
# 
# dat <- merge(dat,ref,by=c("fzone","doy"))
# dat[,`:=`(dog_mean_anom = dog_mean - dog_u)]
# 
# 
# # dat <- merge(dat, 
# #              dflow, 
# #              by=c("time", "year","month","doy"))
# 
# dat %>% 
#   ggplot(aes(time, dog_mean_anom))+
#   # geom_point() + 
#   geom_smooth(method='gam') + 
#   facet_wrap(~fzone, scales='free')
# 
# 
# dat_n <- merge(dat[fzone%in%c("NorthZone_gt10","NorthZone_lt10")], 
#       dflow[site=="torrens_holbrook"], 
#       by=c("time","year","month","doy"))
# n1 <- bam(dog_mean_anom ~ s(dv,fzone,bs='fs',k=5), 
#           data=dat_n,
#           discrete= T,
#           select=T)
# n2 <- bam(dog_mean_anom ~ s(dv_max5d,fzone,bs='fs',k=5), 
#           data=dat_n,
#           discrete= T,
#           select=T)
# n3 <- bam(dog_mean_anom ~ s(dv_max15d,fzone,bs='fs',k=5), 
#           data=dat_n,
#           discrete= T,
#           select=T)
# n4 <- bam(dog_mean_anom ~ s(dv_max30d,fzone,bs='fs',k=5), 
#           data=dat_n,
#           discrete= T,
#           select=T)
# bbmle::BICtab(n1,n2,n3,n4)
# summary(n3)
# summary(n1)
# 
# 
# 
# n4 <- bam(dog_mean_anom ~ s(dv_max30d,fzone,bs='fs',k=5), 
#           data=dat_n,
#           discrete= T,
#           select=T)
# n5 <- bam(dog_mean_anom ~ s(dv_sum30d,fzone,bs='fs',k=5), 
#           data=dat_n,
#           discrete= T,
#           select=T)
# bbmle::BICtab(n4,n5)
# summary(n5)
# 
# 
# plot(n5)
# curve(10*exp(-1*x), 0,10); abline(v = 1,h=5)
# 
# dat_n %>% 
#   ggplot(aes(dv_sum90d, dog_mean))+
#   geom_point(alpha= 0.1)+ 
#   geom_smooth(method='gam')
# 
# dat_n[dog_mean<=0]
# bam(I(log(dog_mean)) ~ 
#       s(dv_sum90d), 
#       # s(log(dv_sum30d),k=5), 
#     data=dat_n[fzone=='CentralZone_lt10'][dog_mean>=0],
#     discrete= T,
#     select=T) %>% plot()
# 
# 
# names(dat)
# m1 <- bam(dog_mean_anom ~ s(dv,fzone,bs='fs',k=5), 
#           data=dat,
#           discrete= T,
#           select=T)
# m2 <- bam(dog_mean_anom ~ s(dv_max5d,fzone,bs='fs',k=5), 
#           data=dat,
#           discrete= T,
#           select=T)
# m3 <- bam(dog_mean_anom ~ s(dv_max15d,fzone,bs='fs',k=5), 
#           data=dat,
#           discrete= T,
#           select=T)
# m4 <- bam(dog_mean_anom ~ s(dv_max30d,fzone,bs='fs',k=5), 
#           data=dat,
#           discrete= T,
#           select=T)
# summary(m4)
# plot(m4,pages=1)
# 
# bbmle::AICtab(m1,m2,m3,m4)
# 
# dat$dog_mean_anom %>% hist
# dat$dv_max30d %>% hist()
# 
# 
# b1 <- bam(dog_mean ~ 
#             s(year,fzone,bs=('fs')) + 
#             s(doy, fzone, bs=c('fs','cc')) + 
#             s(dv_max30d,fzone,bs='fs',k=5), 
#           data=dat[dog_mean > 0],
#           # discrete= T,
#           select=T)
# summary(b1)
# plot(b1,pages=1)
# 
# 
# gratia::smooth_estimates(b1) %>% 
#   ggplot(aes(dv_max30d, .estimate,color=fzone))+
#   geom_line()
# 


d_ts <- d3 %>% filter(is.na(sdog_mean)==F) %>% 
  filter(time >= ymd("2001-01-01")) %>% 
  filter(is.na(fzone)==F)


fn_get_ts <- function(sel_zone){
  tmp <- d_ts[fzone==sel_zone]
  ts_out <- ts(data = tmp$sdog_mean, 
            frequency = 365, 
            start = c(2001,1))
  return(ts_out)
}

c("CentralZone_gt10", "CentralZone_lt10", "NorthZone_gt10", "NorthZone_lt10", 
  "SouthZone_gt10", "SouthZone_lt10")

d_ts %>% 
  ggplot(aes(time, sdog_mean))+
  geom_line()+
  labs(title = "Smoothed and Gapfilled Dogliotti") +
  facet_wrap(~fzone,ncol=1, 
             scales='free') + 
  theme_linedraw()
ggsave(filename = "figures/mcd43a4_smoothed_gapfilled_dogliotti_by_zone.png", 
       width = 20,
       height = 15,
       units='cm',
       scale = 1.5,
       dpi=350)


d_ts %>% 
  ggplot(aes(time, sdog_mean_anom))+
  geom_line()+
  labs(title = "Smoothed and Gapfilled Dogliotti Anomaly") +
  facet_wrap(~fzone,ncol=1, 
             scales='free') + 
  theme_linedraw()
ggsave(filename = "figures/mcd43a4_smoothed_gapfilled_dogliotti-anom_by_zone.png", 
       width = 20,
       height = 15,
       units='cm',
       scale = 1.5,
       dpi=350)



d_ts[fzone=="SouthZone_gt10"]$fzone %>% table

tmp1 <- fn_get_ts("NorthZone_gt10") 
tmp1 %>% 
  decompose() %>% plot()
class(tmp1)
plot.ts(tmp1)

as.data.table(tmp1)$x %>% plot

vec_zones %>% 
  lapply(., fn_plot_decomposition)
  
x <- vec_zones[6]
