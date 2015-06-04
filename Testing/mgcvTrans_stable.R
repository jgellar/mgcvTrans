# Stable Tests - Run this first!

library(devtools)
dev_mode()
load_all()
load_all("../refund/")
library(mgcv)

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
logf <- function(x) log(x0)
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
ggplot(est1.1, aes(x,f)) +
  geom_line(data=filter(est1.1, type=="Truth"), size=2) +
  geom_line(data=filter(est1.1, type!="Truth"), aes(colour=type))
mse1.1 <- sapply(c("tprs", "adaptive", "log transformed"), function(x) {
  mean((est1.1$f[est1.1$type==x] - est1.1$f[est1.1$type=="Truth"])^2)
})




# VDFR
data(sofa)
fit2.1a <- pfr(death ~ lf.vd(SOFA), data=sofa)




est2.1a <- coef(fit2.1a)
