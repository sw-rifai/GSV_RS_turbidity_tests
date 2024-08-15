## SCRATCH ====================
# 2019-08-10
z1 <- d2[
  ,.(val = max(dog,na.rm=T),
     val2 = mean(red,na.rm=T),
     nobs = sum(is.na(dog)==F)),by=.(time,dz)
]
#[,`:=`(time = ymd(paste(year,month,1)))]

z1$val %>% hist
z1$val2 %>% hist
curve(fn_dogliotti2015(x),0.00001,0.16)
z1[,`:=`(year=year(time))][year==2020][dz=="CentralZone_lt10"] %>% 
  .[val2 > 0.05]

z1[nobs>100] %>% 
  ggplot(aes(time, val2,color=dz))+
  geom_line(alpha=0.5)+
  # geom_smooth(method='gam',
  #             aes(weight = nobs),
  #             formula = y~s(x,bs='ad'), 
  #             method.args = list(family=Gamma(link='log'))) + 
  facet_wrap(~dz,ncol = 2)


time(s2)
d2[sample(.N, 10000)]$red %>% hist(1000)
tmp <- d2[time==ymd("2019-08-10")]
tmp$dog %>% summary
tmp$red %>% min(na.rm=T)
fn_dogliotti2015(0.01)
tmp %>% 
  filter(y > -34.9) %>% 
  filter(dog > 0) %>% 
  filter(dog < 400) %>% 
  ggplot(aes(x,y,fill=log10(dog)))+
  geom_raster()+
  scale_fill_viridis_c(option='H',
                       # limits=c(0.01,0.2),
                       oob=scales::squish)+
  coord_sf()
tmp %>% 
  filter(y > -34.9) %>% 
  pull(red) %>% hist(1000)

time(s2)[1:10]
time(s2)[51:76]
tmp2 <- d2[time==ymd("2019-01-30")]
tmp2 %>% 
  filter(y > -34.9) %>% 
  filter(dog > 0) %>% 
  filter(dog < 400) %>% 
  ggplot(aes(x,y,fill=dog))+
  geom_raster()+
  scale_fill_viridis_c(option='H',
                       limits=c(0,15),
                       oob=scales::squish)+
  coord_sf()

