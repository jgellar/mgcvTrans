# Stable Tests - Run this first!

#library(devtools)
library(mgcv)
devtools::load_all()
library(refund)
library(ggplot2)
library(gridExtra)
library(dplyr)

# Set up data
data(DTI)
DTI1 <- DTI[DTI$visit==1 & complete.cases(DTI),]
f <- function(x) x^2

N <- 500
x <- exp(runif(N, 0, 2*pi))
xtrue <- exp(seq(0,2*pi,length=1000))



############
# dt Basis #
############

# Simple log domain transformation - fit 3 ways
dat1.1 <- data.frame(x=x)
dat1.1$Y <- sin(log(dat1.1$x)) + rnorm(N, sd=.5)
mod.dt1 <- gam(Y ~ s(x), data=dat1.1)
mod.dt2 <- gam(Y ~ s(x, bs="ad"), data=dat1.1)
mod.dt3 <- gam(Y ~ s(x, bs="dt", xt=list(tf="log")), data=dat1.1)

est1.1 <- data.frame(type="Truth", x=xtrue, f=sin(log(xtrue)))
est1.1 <- rbind(est1.1, data.frame(type="tprs", x=xtrue,
                                   f=predict(mod.dt1, newdata=data.frame(x=xtrue))))
est1.1 <- rbind(est1.1, data.frame(type="adaptive", x=xtrue,
                                   f=predict(mod.dt2, newdata=data.frame(x=xtrue))))
est1.1 <- rbind(est1.1, data.frame(type="log transformed", x=xtrue,
                                   f=predict(mod.dt3, newdata=data.frame(x=xtrue))))

p1.1 <- ggplot(est1.1, aes(x,f)) +
  geom_line(data=filter(est1.1, type=="Truth"), size=2) +
  geom_line(data=filter(est1.1, type!="Truth"), aes(colour=type))
p1.1
p1.1 + scale_x_log10()

mse1.1 <- sapply(c("tprs", "adaptive", "log transformed"), function(x) {
  mean((est1.1$f[est1.1$type==x] - est1.1$f[est1.1$type=="Truth"])^2)
})
mse1.1


# Quantile transformation when observations are bunched up near 1
dat1.2 <- data.frame(x=(runif(N, 0, 1)^(1/4)))
dat1.2$y <- sin(2*pi*(dat1.2$x)^6) + rnorm(N, sd=.5)
mod.ecdf1 <- gam(y ~ s(x), data=dat1.2)
mod.ecdf2 <- gam(y ~ s(x, bs="ad"), data=dat1.2)
mod.ecdf3 <- gam(y ~ s(x, bs="dt", xt=list(tf="log")), data=dat1.2)
mod.ecdf4 <- gam(y ~ s(x, bs="dt", xt=list(tf="ecdf")), data=dat1.2)

xtrue <- seq(0.001,1,by=.001)
est1.2 <- data.frame(type="Truth", x=xtrue, f=sin(2*pi*xtrue^6))
est1.2 <- rbind(est1.2, data.frame(type="tprs", x=xtrue,
                                   f=predict(mod.ecdf1, newdata=data.frame(x=xtrue))))
est1.2 <- rbind(est1.2, data.frame(type="adaptive", x=xtrue,
                                   f=predict(mod.ecdf2, newdata=data.frame(x=xtrue))))
est1.2 <- rbind(est1.2, data.frame(type="log transformed", x=xtrue,
                                   f=predict(mod.ecdf3, newdata=data.frame(x=xtrue))))
est1.2 <- rbind(est1.2, data.frame(type="ecdf", x=xtrue,
                                   f=predict(mod.ecdf4, newdata=data.frame(x=xtrue))))

p1.2 <- ggplot(est1.2, aes(x,f)) +
  geom_line(data=filter(est1.2, type=="Truth"), size=2) +
  geom_line(data=filter(est1.2, type!="Truth"), aes(colour=type)) +
  scale_y_continuous(limits=c(-1.25,1.25))
p1.2


# Evenly spaced data, but f is more variable at one end than the other
N <- 500
dat3 <- data.frame(x=runif(N))
dat3$Y1 <- sin(2*pi*log(dat3$x)) + rnorm(N, sd=.5)
dat3$Y2 <- sin(2*pi*(dat3$x)^6) + rnorm(N, sd=.5)
xtrue  <- seq(.001, 1, by=.001)
ytrue1 <- sin(2*pi*log(xtrue))
ytrue2 <- sin(2*pi*xtrue^6)

fit3.1a <- gam(Y1 ~ s(x), data=dat3)
fit3.1b <- gam(Y1 ~ s(x, bs="ad"), data=dat3)
fit3.1c <- gam(Y1 ~ s(x, bs="dt", xt=list(tf="log")), data=dat3)
fit3.1d <- gam(Y1 ~ s(x, bs="dt", xt=list(tf="ecdf")), data=dat3)
fit3.1e <- gam(Y1 ~ s(x, bs="dt", xt=list(tf=function(x) x^2)), data=dat3)
fit3.1f <- gam(Y1 ~ s(x, bs="dt", xt=list(tf=function(x) x^6)), data=dat3)
mse3.1 <- sapply(list(fit3.1a, fit3.1b, fit3.1c, fit3.1d, fit3.1e, fit3.1f),
                 function(fit) {
                   est <- predict(fit, newdata=data.frame(x=xtrue), type="terms")
                   mean((est-ytrue1)^2)
                 })
names(mse3.1) <- c("tprs", "adaptive", "log", "ecdf", "x^2", "x^6")

fit3.2a <- gam(Y2 ~ s(x), data=dat3)
fit3.2b <- gam(Y2 ~ s(x, bs="ad"), data=dat3)
fit3.2c <- gam(Y2 ~ s(x, bs="dt", xt=list(tf="log")), data=dat3)
fit3.2d <- gam(Y2 ~ s(x, bs="dt", xt=list(tf="ecdf")), data=dat3)
fit3.2e <- gam(Y2 ~ s(x, bs="dt", xt=list(tf=function(x) x^2)), data=dat3)
fit3.2f <- gam(Y2 ~ s(x, bs="dt", xt=list(tf=function(x) x^6)), data=dat3)
mse3.2 <- sapply(list(fit3.2a, fit3.2b, fit3.2c, fit3.2d, fit3.2e, fit3.2f),
                 function(fit) {
                   est <- predict(fit, newdata=data.frame(x=xtrue), type="terms")
                   mean((est-ytrue2)^2)
                 })
names(mse3.2) <- c("tprs", "adaptive", "log", "ecdf", "x^2", "x^6")




# AR

N <- 1000
d <- tibble(agecat = sample.int(8, N, replace = TRUE)) %>%
  mutate(mu = sin(pi*agecat/7),
         Y = mu + rnorm(N, sd = .1))
fit4.1 <- gam(Y ~ s(agecat, k = 4), data = d)
fit4.2 <- gam(Y ~ s(agecat, bs = "ar"), data = d)





# VDFR
data(sofa)
fit2.1a <- pfr(death ~ lf.vd(SOFA), data=sofa)




est2.1a <- coef(fit2.1a)
