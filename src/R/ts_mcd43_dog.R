library(tidyverse)
library(data.table)
library(mgcv)
library(gratia)

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

fn_kalman <- function(sel_dat){
  y <- sel_dat$dog_mean
  
  ## Set constant parameters:
  dt <- ct <- matrix(0) 
  Zt <- Tt <- matrix(1)
  a0 <- y[1:50] %>% mean(na.rm=T) %>% as.numeric()            # Estimation of the first year flow 
  P0 <- matrix(100)     # Variance of 'a0'
  
  ## Estimate parameters:
  fit.fkf <- optim(c(HHt = var(y, na.rm = TRUE) * .5,
                     GGt = var(y, na.rm = TRUE) * .5),
                   fn = function(par, ...)
                     -fkf(HHt = matrix(par[1]), GGt = matrix(par[2]), ...)$logLik,
                   yt = rbind(y %>% as.numeric), #rbind(y), 
                   a0 = a0, P0 = P0, dt = dt, ct = ct,
                   Zt = Zt, Tt = Tt)
  
  
  
  ## Filter Nile data with estimated parameters:
  fkf.obj <- fkf(a0, P0, dt, ct, Tt, Zt, HHt = matrix(fit.fkf$par[1]),
                 GGt = matrix(fit.fkf$par[2]), yt = rbind(y %>% as.numeric)
  )
  sel_dat$kdog <- fkf.obj$att[1,]
  return(sel_dat)
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

x <- vec_zones[2]
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
    mutate(fzone = x) %>% 
    mutate(dog_gf = coalesce(dog_mean, sdog_mean))

  
  tmp_dat <- fn_kalman(tmp_dat)
  tmp_dat[tmp_dat$kdog < 0]$kdog <- 0.001
  tmp_dat <- tmp_dat %>% 
    mutate(kdog_gf = coalesce(dog_mean,kdog))
  return(tmp_dat)
}


d2 <- vec_zones %>% lapply(., fn_smooth)
d2 <- rbindlist(d2) %>% 
  mutate(fzone = factor(fzone))
d2[,`:=`(ddate = decimal_date(time))]
d2 <- d2[year >= 2002]


# calculate anomaly ===========================================

m0 <- bam(kdog_gf ~ 
  #log(dog_gf) ~ 
            s(doy,fzone,bs=c("fs","cc"),k=6),
          # s(doy,by=fzone,bs='cc'), 
          data=d2,
          weights = d2$count,
          family = Gamma(link='log'),
          select = T, 
          discrete = T)
summary(m0)
# plot(m0)
# draw(m0)

gratia::smooth_estimates(m0) %>% 
  ggplot(aes(doy, .estimate,color=fzone))+
  geom_line()

ref <- expand.grid(doy = 1:365,
                   fzone = levels(dat$fzone)) %>% setDT()
ref <- ref %>% mutate(dog_u = predict(m0, newdata=., type='response'))

d3 <- merge(d2,ref,by=c("fzone","doy"))
d3[,`:=`(kdog_gf_anom = kdog_gf - dog_u)]

library(ggfortify)
vec_years <- 2002:2022
fn_acf <- function(sel_year){
  d3[fzone == "CentralZone_lt10"][year(time)==sel_year][order(time)]$kdog_gf_anom %>% 
    acf(plot = F, 
        ci.type='ma') %>% 
    fortify() %>% 
    mutate(year=sel_year)
  
}
d_acf <- vec_years %>% 
  lapply(., fn_acf) %>% 
  rbindlist()



d_acf %>% 
  ggplot(aes(Lag,ACF))+
  geom_line(data=d_acf[,.(ACF_u = mean(ACF)),by=Lag],
            aes(Lag, ACF_u),
            inherit.aes = F, 
            color = 'grey30')+
  geom_line(color='#cf0000')+
  labs(x = "Lag days",
       y = "Autocorrelation", 
       title = "Dogliotti Turbidity, Central Zone, < 10 m") +
  facet_wrap(~year) +
  theme_linedraw()+
  theme(panel.grid = element_blank())
ggsave(filename = "figures/ACF_kdog_anom_by_year_CentralZone_lt10.png",
       width = 22,
       height =15,
       units = 'cm')


vv <- d3[fzone == "CentralZone_lt10"][year <= 2022] %>% 
  group_by(year) %>% 
  filter(time >= time[kdog_gf == max(kdog_gf)])
vv <- vv %>% ungroup() %>% as.data.table()
vv %>% 
  mutate(doy = yday(time)) %>% 
  ggplot(aes(doy, kdog_gf_anom))+
    geom_line()+
    geom_line(aes(doy, dog_mean-dog_u),col='red') +
    labs(y = "Dogliotti Anomaly",
         x = "Day of year",
         title = "Central < 10m, post-annual-max",
         subtitle = "black: Kalman filter gapfill") +
    theme_classic() + 
    facet_wrap(~year,scales='free')
ggsave(filename = "figures/timeseries-post-max_kdog_anom_by_year_CentralZone_lt10.png",
       width = 22,
       height =15,
       units = 'cm')


names(d3)

d3[fzone == "CentralZone_lt10"][year == 2021] %>% 
  ggplot(aes(doy, kdog_gf_anom))+
  geom_point(aes(doy,dog_mean - dog_u,color=count),
            alpha =1) + 
  geom_line() + 
  scale_color_viridis_c(option='H',
                        direction = -1) + 
  labs(x = "day of year",
       y = "Dogliotti Anomaly",
       color = "Pixel Count",
       title = "2021, Central Zone, < 10 m", 
       subtitle = "black: Kalman filter gapfill, red: MCD43 zonal mean") + 
  theme_linedraw()
ggsave(filename = "figures/timeseries_kdog_anom_2022_CentralZone_lt10.png",
       width = 22,
       height =15,
       units = 'cm')



vv <- d3[fzone == "CentralZone_lt10"][year <= 2022] %>% 
  group_by(year) %>% 
  filter(time >= time[kdog_gf == max(kdog_gf)])
vec_years <- 2002:2022

fn_acf <- function(sel_year){
  vv[year==sel_year][order(time)]$kdog_gf_anom %>% 
    acf(plot = F, 
        ci.type='ma', 
        lag.max = 60) %>% 
    fortify() %>% 
    mutate(year=sel_year)
  
}
d_acf2 <- vec_years %>% 
  lapply(., fn_acf) %>% 
  rbindlist()
d_acf2 %>% 
  ggplot(aes(Lag,ACF))+
  geom_ribbon(aes(Lag,ymax=upper,ymin=lower),
              alpha=0.2)+ 
  geom_line(col="#cf0000")+
  geom_line(data=d_acf2[,.(ACF_u = mean(ACF)),by=Lag],
            aes(Lag, ACF_u),
            inherit.aes = F, 
            color = 'grey30')+
  labs(x = "Lag Days",
       y = "Dogliotti Anomaly Autocorrelation",
       title = "Central Zone, < 10 m",
       subtitle = "black: average, red: target year") +
  facet_wrap(~year) +
  theme_linedraw()
ggsave(filename = "figures/ACF-post-max_kdog_anom_by_year_CentralZone_lt10.png",
       width = 22,
       height =15,
       units = 'cm')



dat_s2 <- arrow::read_parquet("outputs/S2-Level2A_AMC_red.parquet") %>% 
  mutate(year = year(time),
         month = month(time), 
         doy = yday(time)) %>% 
  mutate(fzone = factor(dz)) %>% 
  mutate(dog_mean = fn_dogliotti(mean),
         nechad_mean = fn_nechad(mean))

dat_s2[fzone=='CentralZone_lt10'][count >= 10000] %>% 
  ggplot(aes(time, dog_mean,color = count))+
  geom_point() + 
  geom_smooth(span = 0.05) + 
  geom_smooth(method = 'gam',
              formula = y~s(x,bs='tp',k=36),
              col = "#cf0000") + 
  scale_color_viridis_c(option="H",direction = -1)


d3[fzone=='CentralZone_lt10'][kdog <= 100][count >= 10] %>% 
  ggplot(aes(time, kdog_gf_anom,color = count))+
  geom_point() + 
  # geom_line(aes(time, dog_u),col='grey30',lwd=2) + 
  # geom_smooth(span = 0.05) + 
  # geom_smooth(method = 'gam',
  #             formula = y~s(x,bs='tp',k=36),
  #             col = "#cf0000") + 
  scale_color_viridis_c(option="H",direction = -1)

d3[fzone=='CentralZone_lt10'][kdog <= 100][count >= 10][time>=ymd("2019-01-01")]

rbindlist(list(
  d3[fzone=='CentralZone_lt10'][kdog <= 100][count >= 10][time>=ymd("2019-01-01")] %>% 
  select(time, kdog_gf,count) %>% 
  rename(dog = kdog_gf) %>% 
  mutate(product = "MCD43"),
  dat_s2[fzone=='CentralZone_lt10'][count >= 10000] %>% 
  select(time,dog_mean,count) %>% 
  rename(dog = dog_mean) %>% 
  mutate(product = "S2_L2A"))) %>% 
  ggplot(aes(time, dog,
             color=product))+
  geom_point() +
  geom_smooth(aes(weight = count),
              method='gam',
              formula = y~s(x,bs='tp',k=50), 
              lwd = 1.25,
              alpha = 0.3) + 
  geom_smooth(method='lm', aes(weight = count),
              se = F,
              color = 'black',
              lwd =2.5) + 
  geom_smooth(method='lm', aes(weight = count),
              se = F, 
              lty =1) + 
  scale_color_brewer(type = "qual", palette = 2)+
  theme_linedraw()
ggsave(filename = "figures/compare_mcd43_s2_dog_CentralZone_lt10.png",
       width = 22,
       height =15,
       units = 'cm',
       dpi = 350)



dev.new()
d3[fzone=='CentralZone_lt10'][kdog <= 100][count >= 10][
  time>=ymd("2021-06-01")
  ][time <= ymd('2022-06-01')] %>% 
  ggplot(aes(time, q95))+
  geom_point()



 d3[fzone == "CentralZone_lt10"][year(time)==2007][order(time)] %>% 
  ggplot(aes(time,dog_gf_anom))+
  geom_line()

d3[fzone == "CentralZone_lt10"][order(time)] %>%
  mutate(doy = yday(time)) %>% 
  ggplot(aes(doy,dog_mean))+
  geom_line(aes(doy,dog_gf),col='red') + 
  geom_line()+
  facet_wrap(~year, 
             scales = 'free')


tmp1 <- d2[fzone == "CentralZone_lt10"][order(time)]
tmp1$dog_mean %>% is.na %>% table


m0 <- gam(dog_mean ~ s(ddate),
           data=tmp1[ddate >=2002][1:1000,],
          select = T)
plot(m0)

m1 <- gamm(dog_mean ~ s(ddate),
     correlation = corAR1(), 
     select = T,
     method = "REML",
     data=tmp1[ddate >=2010][1:3000,])
plot(m1$gam)
m1$lme

nlme(dog_mean ~ 1,
     correlation = corAR1(),
  data=tmp1[ddate >=2002][1:1000,])

as_tsibble(tmp1, index = time) %>% 
  select(sdog_mean) %>% 
  PACF(sdog_mean, 
      # type='covariance') %>% 
  autoplot()

data("nycflights")

vec_d <- tmp1$dog_mean
tmp1$time

tmp1$time %>% yday %>% first
ts_dog <- zoo::zooreg(vec_d, 
                      start = 87,
                      order.by=tmp1$time, 
                      frequency = 365,
                      deltat = 1/365)
plot(ts_dog)
ts_dog
?na.approx
ts_dog2 <- na.fill(ts_dog,'extend')
decompose(ts_dog2)
pacf(ts_dog2)


tmp1_p <- arima(ts_dog)
predict(tmp1_p)$pred %>% plot

summary(tmp1$dog_mean)
tmp1_p <- arima(tmp1$dog_mean) %>% predict
tmp1_p$pred %>% plot
