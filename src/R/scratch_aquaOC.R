## SCRATCH ====================
ref <- d2[time <= ymd('2019-06-01')][
  ,.(val_max = max(dog,na.rm=T),
     val_med = median(dog,na.rm=T),
     val_min = min(dog,na.rm=T),
     nobs = sum(is.na(dog)==F)),by=.(year,month,dz)
][,.(val_max = median(val_max,na.rm=T),
     val_med = median(val_med,na.rm=T),
     val_min = min(val_min,na.rm=T)),by=.(month,dz)]

# 2019-08-10
z1 <- d2[
  ,.(val_max = max(dog,na.rm=T),
     val_med = median(dog,na.rm=T),
     val_min = min(dog,na.rm=T),
     nobs = sum(is.na(dog)==F)),by=.(year,month,dz)
][,`:=`(time=ymd(paste(year,month,1)))]

z1[dz %in% levels(z1$dz)[1:4]] %>% 
  merge(., ref, by=c("month","dz"),
        suffix = c("","_ref")) %>% 
  # mutate(val = val+0.001) %>% 
  # filter(time >= ymd("2015-01-01")) %>% 
  # filter(time <= ymd("2021-01-01")) %>% 
  ggplot(aes(time, val_med))+
  geom_ribbon(aes(ymin=val_min_ref,ymax=val_max_ref),
              alpha=0.8,
              fill='black',
              lty=0)+
  geom_ribbon(aes(ymin=val_min,ymax=val_max),
              alpha=0.5,
              fill='red')+
  geom_line(aes(time,val_med_ref))+
  geom_line(col='red')+
  # geom_smooth(method='gam',
  #             aes(weight = nobs), 
  #             formula = y~s(x,bs='ad',k=30),
  #             method.args = list(family=Gamma(link='log'))) +
  scale_color_viridis_d(option='H') + 
  scale_x_date(date_breaks = '1 year',
               date_labels = '%Y') +
  facet_wrap(~dz,ncol = 2)

z1[dz %in% levels(z1$dz)[1:4]] %>% 
  merge(., ref, by=c("month","dz"),
        suffix = c("","_ref")) %>% 
  mutate(anom_med = val_max - val_max_ref) %>% 
  group_by(year,dz) %>% 
  summarize(a_med = sum(anom_med)) %>% 
  ungroup() %>% 
  ggplot(aes(year,a_med,color=dz))+geom_line()+
  scale_color_brewer(type='qual',palette = 2)+
  scale_x_continuous(breaks = 2002:2024)




zones %>% ggplot(aes(x,y,fill=dz))+geom_tile()+coord_sf()+
  scale_fill_viridis_d(option='H')
