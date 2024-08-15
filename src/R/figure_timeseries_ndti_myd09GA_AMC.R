pacman::p_load(lubridate, terra, tidyverse, data.table,
               mgcv)

dat <- fread("data/gee_SAWater_turbidity/MYD09GA_NDTI_medianCountSD_2003-01-01_2021-12-31.csv") %>% 
  rename(ndti = median) %>% 
  .[order(zone,date)] %>% 
  select(-`.geo`) %>% 
  mutate(year=year(date),
         month=month(date), 
         fzone = factor(zone, ordered = T,
         levels = c("North Zone","Central Zone","South Zone"),    
         labels = c("North Zone","Central Zone","South Zone"))) %>% 
  .[]
dat[,`:=`(ndti_4w = frollmean(ndti, n=4, fill=NA, align='right'))]
dat[,`:=`(ndti_12w = frollmean(ndti, n=12, fill=NA, align='right'))]
dat[,`:=`(ndti_52w = frollmean(ndti, n=52, fill=NA, align='right'))]

dat[date %between% c(ymd("2009-01-01"),ymd("2010-12-31"))] %>% 
  ggplot(aes(date, ndti, color=count))+
  geom_hline(aes(yintercept = 0),
    color = 'grey70')+ 
  geom_ribbon(aes(ymin = ndti - stdDev, 
                    ymax = ndti + stdDev),
              alpha = 0.25, 
              color = NA) + 
  geom_point() + 
  scale_color_viridis_c() +
  labs(x = NA, y="NDTI",
    color = "Pixel Obs") + 
  facet_wrap(~zone, ncol = 1)+
  theme_linedraw()+
  theme(panel.grid = element_blank())
  
dat[date %between% c(ymd("2009-01-01"),ymd("2010-12-31"))] %>% 
  ggplot(aes(count, ndti))+
  geom_point()+
facet_wrap(~zone, ncol = 1)

dat %>% 
ggplot(aes(date, ndti_4w, color=zone))+
  geom_line()


dat %>% 
  ggplot(aes(date, ndti_12w, color=zone))+
  geom_line()

dat %>% 
  ggplot(aes(date, ndti_52w, color=zone))+
  geom_line()


dat %>% 
  select(ndti,zone,date) %>% 
  dcast(data=., date + zone ~ ndti)

dat %>% 
  select(ndti, zone, date) %>% 
  pivot_wider(values_from = ndti, 
              names_from = c(zone,date))



library(mgcViz)
## start with seasonal component and zone
m1 <- gam(ndti ~ s(month, by = zone, 
                   bs=c('cc','fs')), 
          data=dat, 
          select = T)
plot(m1)

getViz(m1) %>% 
  plot() %>% 
  print(pages=1)

m2 <- gam(ndti ~ zone + 
                   s(month, 
                   bs=c('cc')), 
          data=dat, 
          select = T)

bbmle::AICctab(m1,m2)
## Lower AIC on m2 suggests that each zone does not need its own monthly smoooth.


## Add year component0
m3 <- gam(ndti ~ year + 
            zone +  
                   s(month,  
                   bs=c('cc')), 
          data=dat, 
          select = T)
bbmle::AICctab(m2, m3)
summary(m3)
getViz(m3) %>% plot(allTerms = T) %>% print(pages = 1)
## Relative to zone differences, the linear year effect is extremely weak

## Examine year as a nonlinear effect
m4 <- gam(ndti ~ 
            s(year) + 
            zone +
            s(month,
            bs=c('cc')), 
          data=dat,
          select = T)
getViz(m4) %>% plot(allTerms = T) %>% print(pages = 1)
summary(m4)
bbmle::AICctab(m3, m4)
## year is better modeled as a nonlinear effect, than a linear effect


## Examine nonlinear year effect by site
m5 <- gam(ndti ~ 
            s(year, zone, bs='fs') + 
            s(month,
              bs=c('cc')), 
          data=dat,
          select = T)
getViz(m5) %>% plot(allTerms = T) %>% print(pages = 1)
summary(m5)
bbmle::AICctab(m4, m5)
## nonlinear year by site is better, but the AIC difference is very small
## so we will drop the year x zone interaction


## Examine number of knots on year
kn <- dat$year %>% unique %>% length
m6 <- gam(ndti ~ 
            s(year, k = kn, bs = 'ts') + 
            zone +
            s(month,
              bs=c('cc')), 
          data=dat,
          select = T)
getViz(m6) %>% plot(select = 1) %>% print(pages = 1)
summary(m6)
bbmle::AICctab(m5, m6)
## 
