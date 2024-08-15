library(ggfortify)
library(tidyverse)
library(data.table)
library(mgcv)
library(gratia)
library(cols4all)
library(patchwork)
library(FKF)

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

fn_smooth <- function(x, sel_data){
  tmp_dat <- sel_data[fzone==x]
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

fn_theme <- function (base_size = 14, base_family = "", base_line_size = base_size/22, 
                      base_rect_size = base_size/22) {
  half_line <- base_size/2
  theme_bw(base_size = base_size, base_family = base_family, 
           base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace% 
    theme(axis.text = element_text(colour = "black", size = rel(0.8)), 
          axis.ticks = element_line(colour = "black", 
                                    linewidth = rel(0.5)),
          axis.ticks.length = unit(-0.1,'cm'),
          panel.border = element_rect(fill = NA, colour = "black", 
                                      linewidth = rel(1)), panel.grid = element_line(colour = "black"), 
          panel.grid.major = element_blank(), #element_line(linewidth = rel(0.1)), 
          panel.grid.minor = element_blank(), #element_line(linewidth = rel(0.05)), 
          strip.background = element_rect(fill = "transparent",
                                          color = 'transparent'), 
          strip.text = element_text(colour = "black", size = rel(0.8), 
                                    margin = margin(0.8 * half_line, 0.8 * half_line, 
                                                    0.8 * half_line, 0.8 * half_line)), 
          legend.position = c(0.5,0.025),
          legend.justification = c(0.5,0.025),
          legend.text = element_text(size = rel(0.8)),
          legend.direction = 'horizontal',
          legend.background =  element_rect(fill='transparent', 
                                            color='transparent'),
          complete = TRUE) 
}


# Import Data  ==============================================
dat_s2 <- arrow::read_parquet("outputs/S2-Level2A_AMC_red.parquet") %>% 
  mutate(year = year(time),
         month = month(time), 
         doy = yday(time)) %>% 
  mutate(fzone = factor(dz)) %>% 
  mutate(dog_mean = fn_dogliotti(mean),
         nechad_mean = fn_nechad(mean))

dat_m <- arrow::read_parquet("outputs/MCD43A4_AMC_red.parquet") %>% 
  mutate(year = year(time),
         month = month(time), 
         doy = yday(time)) %>% 
  mutate(ddate = decimal_date(time)) %>% 
  mutate(fzone = factor(dz)) %>% 
  mutate(dog_mean = fn_dogliotti(mean),
         nechad_mean = fn_nechad(mean))

# Smooth time series ==========================================
vec_zones <- dat_s2$fzone %>% levels()
d_time <- data.table(time = seq(dat_s2[mean>0]$time %>% min,
                                to = dat_s2[mean>0]$time %>% max,
                                by='1 day'))


d2 <- vec_zones %>% lapply(., fn_smooth, sel_data = dat_s2)
d2 <- rbindlist(d2) %>% 
  mutate(fzone = factor(fzone))
d2[,`:=`(ddate = decimal_date(time))]
d2 <- d2[year >= 2002]

# calculate anomaly ===========================================

m0 <- bam(kdog_gf ~ 
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
                   fzone = levels(dat_s2$fzone)) %>% setDT()
ref <- ref %>% mutate(dog_u = predict(m0, newdata=., type='response'))

d3 <- merge(d2,ref,by=c("fzone","doy"))
d3[,`:=`(kdog_gf_anom = kdog_gf - dog_u)]
d3[,`:=`(count_max = max(count,na.rm=T)),by=fzone]
d3[,`:=`(yr_q = quarter(time))]



fn_acf <- function(x){
  x[order(time)]$kdog_gf_anom %>% 
    acf(plot = F, lag.max = 30,
        ci.type='ma') %>% 
    fortify()
}
d_acf <- d3[count > (count_max*0.1)][,fn_acf(.SD),by=.(fzone,year,yr_q)]

d_acf %>% 
  ggplot(aes(Lag, ACF,color=factor(yr_q)))+
  geom_hline(aes(yintercept = 0),lwd = 1) + 
  # geom_ribbon(aes(ymin=lower,ymax=upper))+
  geom_line() + 
  facet_grid(fzone ~ year)


names(d3)
d3
d3[fzone=="CentralZone_lt10"][year%in%c(2019:2021)] %>% 
  ggplot(aes(doy,fn_dogliotti(q95),color=factor(year)))+
  geom_point()+
  geom_smooth(se=F,
              aes(weight = count), 
              method='gam',
              formula = y~s(x,bs='ad',k=30))


library(cols4all)
p1 <- d3 %>% 
  ggplot(aes(time, mean,color=fzone))+
  geom_smooth(method='lm',se=F)+
  scale_color_discrete_c4a_cat(palette = "viridis.turbo") + 
  labs(title = "S2-L2A") + 
  coord_cartesian(ylim = c(0,0.05))


p2 <- dat_m[time %between% (d3$time %>% range())] %>% 
  ggplot(aes(time, mean,color=fzone))+
  geom_smooth(method='lm',se=F)+
  scale_color_discrete_c4a_cat(palette = "viridis.turbo") +
  labs(title= "MCD43")+
  coord_cartesian(ylim=c(0,0.05))


d3$ddate
d3$dog_mean
m_s2 <- gam(dog_mean ~ 
              fzone + 
              
              s(ddate, fzone, bs=c("fs","ad")) + 
              s(doy,fzone, bs= c("fs", 'cc')), 
            data=d3,
            weights = count,
            select=F, 
            method='REML')
summary(m_s2)

draw(m_s2)
gratia::smooth_estimates(m_s2) %>%
 filter(.smooth == "s(ddate,fzone)") %>% 
  ggplot(aes(ddate, estimate, color= ))



m_m <- gam(dog_mean ~ s(ddate, fzone, bs='fs') + 
              s(doy,fzone, bs= c("fs", 'cc')), 
            data=dat_m,
            weights = count,
            select=T, 
            method='REML')
draw(m_m)


m_2 <- bam(dog_mean ~ 
             fzone + 
             s(ddate, fzone, bs=c("fs","ad"),k=30) + 
             s(doy,fzone, bs= c("fs", 'cc'),k=12), 
           data=dat_m[year >= 2003],
           weights = count,
           select=T, 
           discrete = T)
draw(m_2)



# Weighted MCD43 anomaly =========================
m_ref_mcd <- bam(dog_mean ~ 
            fzone + 
            s(doy,fzone,bs=c("fs","cc"),k=6),
          data=dat_m[dog_mean > 0],
          weights = count,
          family = Gamma(link='log'),
          select = T, 
          discrete = T)
summary(m_ref_mcd)
draw(m_ref_mcd)

dat_m <- dat_m %>% 
  mutate(dog_u = predict(m_ref_mcd,type='response',newdata=.)) %>% 
  mutate(dog_anom = dog_mean - dog_u)

dat_m2 <- dat_m[is.na(dog_mean)==F][,.(
            val = weighted.mean(dog_anom, count),
             val_sd = sd(dog_anom)),
                          by=.(year,month,fzone)][
      ,`:=`(d_ym = ymd(paste(year,month,1)))
                          ][]
dat_m2 <- dat_m2[order(fzone,d_ym)][,`:=`(val12 = frollmean(val,n=12,align='right')), 
                          by = fzone]

dat_m2 %>% 
  filter(d_ym >= ymd("2002-12-01")) %>% 
  ggplot(aes(d_ym, val12,color=fzone))+
  geom_hline(aes(yintercept = 0),color='black') +
  geom_line(color='black',
            aes(group=fzone),
            lwd = 1.5) + 
  geom_line() + 
  scale_color_brewer(type='qual',palette = 2) + 
  labs(x = NULL, 
       y = "12 mo mean Dog. FNU anom.", 
       color = NULL) + 
  scale_x_date(date_breaks = "2 years", 
               date_labels = "%Y", 
               expand = c(0,0)) + 
  fn_theme() + 
  theme(legend.position = 'bottom')
ggsave(filename = paste0('figures/mcd43a4_weightedDog12mo_timeseries_',Sys.Date(),'.png'),
       width = 30,
       height = 15,
       scale = 0.75,
       units = 'cm',
       dpi = 350)
