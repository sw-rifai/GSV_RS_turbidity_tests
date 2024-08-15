pacman::p_load(lubridate, terra, tidyverse, data.table,
               mgcv)

dat <- fread("data/gee_SAWater_turbidity/MYD09GA_NDTI_medianCountSD_2003-01-01_2023-07-01.csv") %>% 
  rename(ndti = median) %>% 
  .[order(zone,date)] %>% 
  select(-`.geo`) %>% 
  mutate(year=year(date),
         month=month(date), 
         fzone = factor(zone, ordered = T,
                        levels = c("North Zone","Central Zone","South Zone"),    
                        labels = c("North Zone","Central Zone","South Zone"))) %>% 
  .[]
dat[,ndti_4w := frollmean(ndti, n=4, fill=NA, align='right')]
dat[,ndti_12w := frollmean(ndti, n=12, fill=NA, align='right')]
dat[,ndti_52w := frollmean(ndti, n=52, fill=NA, align='right')]

dat2 <- fread("data/gee_SAWater_turbidity/MOD09GA_NDTI_medianCountSD_2000-03-01_2023-07-01.csv") %>% 
  rename(ndti = median) %>% 
  .[order(zone,date)] %>% 
  select(-`.geo`) %>% 
  mutate(year=year(date),
         month=month(date), 
         fzone = factor(zone, ordered = T,
                        levels = c("North Zone","Central Zone","South Zone"),    
                        labels = c("North Zone","Central Zone","South Zone"))) %>% 
  .[]
dat2[,ndti_4w := frollmean(ndti, n=4, fill=NA, align='right')]
dat2[,ndti_12w := frollmean(ndti, n=12, fill=NA, align='right')]
dat2[,ndti_52w := frollmean(ndti, n=52, fill=NA, align='right')]


dat <- rbind(dat,dat2)


dat[date %between% c(ymd("2003-01-01"),ymd("2004-01-01"))][zone=="South Zone"] %>% 
  ggplot(aes(date, ndti_4w,color=product,size = count))+
  geom_point()



dat[date %between% c(ymd("2003-01-01"),ymd("2004-01-01"))] %>% 
  # .[zone=="North Zone"] %>% 
  dcast(data=., zone + date  ~ product, 
        value.var = "ndti_12w",
        fun.aggregate = median,
        drop = T) %>% 
  set_names(tolower(names(.))) %>% 
  ggplot(aes(myd09ga,mod09ga))+
  geom_point() + 
  geom_smooth(method='lm') + 
  geom_abline(col='#cf0000') + 
  coord_equal() + 
  facet_wrap(~zone)
   
dat[date %between% c(ymd("2003-01-01"),ymd("2004-01-01"))][zone=="South Zone"] %>% 
  pivot_wider(values_from = c("ndti","date"),
              names_from = c("product")) %>% 
  select(ndti_MYD09GA, ndti_MOD09GA)
  # select(MYD09GA, MOD09GA) %>% drop_na()
  ggplot(aes(ndti_MYD09GA, ndti_MOD09GA))+
  geom_point()


ChickWeight = as.data.table(ChickWeight)
setnames(ChickWeight, tolower(names(ChickWeight)))
DT <- melt(as.data.table(ChickWeight), id=2:4) # calls melt.data.table

names(DT)
# dcast is an S3 method in data.table from v1.9.6
dcast(DT, time ~ variable, fun=mean) # using partial matching of argument
dcast(DT, diet ~ variable, fun=mean)
dcast(DT, diet+chick ~ time, drop=FALSE)
dcast(DT, diet+chick ~ time, drop=FALSE, fill=0)