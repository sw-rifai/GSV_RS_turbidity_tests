pacman::p_load(lubridate, terra, tidyverse, data.table,
               mgcv, mgcViz,sf)

kd <- fread("data/gee_SAWater_turbidity/MOD09GA_Kd490_medianCountSD_daily_2000-03-01_2023-07-01.csv") %>% 
 rename(kd = median) %>% 
  .[order(zone,date)] %>% 
  select(-`.geo`) %>% 
  mutate(year=year(date),
         month=month(date), 
         fzone = factor(zone, ordered = T,
                        levels = c("North Zone","Central Zone","South Zone"),    
                        labels = c("North Zone","Central Zone","South Zone"))) %>% 
  .[]
kd[,kd_4w := frollmean(kd, n=4, fill=NA, align='right')]
kd[,kd_12w := frollmean(kd, n=12, fill=NA, align='right')]
kd[,kd_52w := frollmean(kd, n=52, fill=NA, align='right')]

flow <- fread("data/daily_discharge_torrens_river_seaview_road_bridge.csv")
names(flow) <- c("date","time","dv")
flow[,`:=`(date = dmy(date))]
flow[,`:=`(year=year(date),month=month(date))]

flow[year >= 2016 & year<=2018] %>% 
  ggplot(aes(date, dv))+
  geom_line()+
  geom_point()

flow[dv == max(dv)]

kd[year >= 2015 & year<=2017] %>% 
  # .[zone=="Central Zone"] %>% 
  # [month %between% c(8,11)] %>% 
  ggplot(aes(date, kd, color = zone))+
  # geom_point()+
  geom_smooth(method = 'gam', 
              formula = y~s(x, bs='ad'),
              aes(weight = count))



dat[,`:=`(year=year(date),month=month(date))]
dat[year >= 2016 & year<=2017][zone=="Central Zone"] %>% 
  # [month %between% c(8,11)] %>% 
  ggplot(aes(date, ndti_4w, color = product))+
  geom_point()+
  geom_smooth(method = 'gam', 
              formula = y~s(x, bs='ad'),
              aes(weight = count))


ref <- dat[,.(ndti_u = mean(ndti,na.rm=T)), 
           by=.(month,product,zone)]
ref
ref %>% 
  ggplot(aes(month,ndti_u,color=paste(zone,product)))+
  geom_line()

dat <- merge(dat, ref, by=c("month","product","zone"))

dat[year >= 2016 & year<=2018][zone=="Central Zone"] %>% 
  # [month %between% c(8,11)] %>% 
  ggplot(aes(date, ndti_4w - ndti_u, color = product))+
  geom_point()+
  geom_smooth(method = 'gam', 
              formula = y~s(x, bs='ad'),
              aes(weight = count))
