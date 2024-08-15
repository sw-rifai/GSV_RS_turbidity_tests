library(tidyverse); 
library(gramEvol)

fn_kd490_oli <- function(bg){
  # https://oceancolor.gsfc.nasa.gov/resources/atbd/kd/
  vec_params <- c(-0.9054,	-1.5245,	2.2392,	-2.4777,	-1.1099)
  
  a0 <- vec_params[1]
  a1 <- vec_params[2]
  a2 <- vec_params[3]
  a3 <- vec_params[4]
  a4 <- vec_params[5]
  log10_kd490 <- a0 + a1*log10(bg) + a2*(log10(bg))**2 + a3*(log10(bg))**3 + a4*(log10(bg))**4
  kd490 <- 0.0166 + 10**log10_kd490
  return(kd490)
}

x <- seq(0.25, 2.5,length.out =100)
y <- fn_kd490_oli(x)
y2 <- log(y)
plot(y~x)

library(mgcv)
l1 <- gam(y2~poly(x, 4))
gam(y2 ~ poly(x,4,raw = T)) %>% coef %>% unname() %>% dput
summary(l1)
plot(l1, all.terms=T)

fn <- function(x) predict(l1, newdata=data.frame(x=x)) %>% exp()
fn(5)
curve(fn(x), 0.25,2.5)
curve(fn_kd490_oli(x), col='red',add=T)

plot(y~predict(l1) %>% exp); abline(0,1,col='red')
plot


ruleDef <- list(expr = grule(op(expr, expr), func(expr), var),
                func = grule(log, sqrt),
                op = grule('+', '-','*'),
                var = grule(x))

grammarDef <- CreateGrammar(ruleDef)
grammarDef

set.seed(321)
GrammarRandomExpression(grammarDef, 6)



# cost function
SymRegFitFunc <- function(expr) {
  result <- eval(expr)
  if (any(is.nan(result)))
    return(Inf)
  out <- (mean(log(1 + abs(y - result))))
  # return (mean(log(1 + abs(y - result))))
  # out <- sqrt( mean( (y-result)**2) )
  # out <- abs(mean(y-result))
  return(out)
}



set.seed(321)

library("parallel")
options(mc.cores = 10)
ge <- GrammaticalEvolution(grammarDef,
                           SymRegFitFunc, 
                           terminationCost = 0.1, 
                           # iterations = 100, 
                           max.depth = 5,
                           # optimizer = 'ga', 
                           monitorFunc = print,
                           plapply = mclapply)
ge

plot(y ~ x)
points(eval(ge$best$expressions) ~ x, col = "red", type = "l")
plot(eval(ge$best$expressions) %>% exp~x,type='l')

fn <- function(x) eval(ge$best$expressions)
fn(1)
curve(fn(x),0.25,2.5)

SymRegFitFunc <- function(expr) {
  result <- eval(expr)
  if (any(is.nan(result)))
    return(Inf)
  return (mean(log(1 + abs(y - result))))
}

set.seed(314)
ge <- GrammaticalEvolution(grammarDef, SymRegFitFunc, terminationCost = 0.1, iterations = 2500, max.depth = 5)
ge
## Grammatical Evolution Search Results:
##   No. Generations:  2149 
##   Best Expression:  sin(x) + cos(x + x) 
##   Best Cost:        0.0923240003917875

plot(y)
points(eval(ge$best$expressions), col = "red", type = "l")