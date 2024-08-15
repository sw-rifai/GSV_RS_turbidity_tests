library(tidyverse); 

curve(3.4725 - 6.2114*x  + 2.7792*x**2,
      0.1,1)

x <- runif(100,0.1,1)
x2 <- runif(100,1,5)
y <- 3.4725 - 6.2114*x  + 2.7792*x**2 + rnorm(100,0,0.1)
y <- c(y, rep(0,100))
x3 <- c(x,x2)

plot(y~x3)

lm(y~log(x3)+x3+I(x3**2)) %>% summary()


curve(3.4725 - 6.2114*x  + 2.7792*x**2,
      0.1,5, ylab="kd490", xlab="Blue/Green",ylim=c(0,10))
curve(-1.01908 -1.88839*log(x) + 1.33892*x -0.10601*x**2, 0.1,5, 
      add=T, col='blue')
legend("topright",
       c("Adi 2013","Modified exp decay polynomial"), 
       lwd=c(1,1),col=c("black","blue"))

points(y~x3)
