# Testing....

library(devtools)
dev_mode()
load_all()
load_all("../FDA - Survival/pcox/")

data(sofa)
tmp <- lf(sofa$SOFA)
cl <- tmp$call
cl[[4]] <- "dt"
names(cl)[4] <- "bs"

spec <- eval(cl, envir = tmp$data)
sm <- smoothCon(spec, data = tmp$data, knots = NULL, absorb.cons = T)

data(DTI)
DTI1 <- DTI[DTI$visit==1 & complete.cases(DTI),]
f <- function(x) log(x)
fit.lf <- pfr(pasat ~ lf(cca, k=30, bs="dt", xt=list(tf=f)), data=DTI1)




chk <- te(var1, var2, xt=list(list(tf="s-t", bs="tp"), list(bs="cr")))
chk <- te(var1, var2, xt=list(bs=c("tp","cr")))


# dt testing - Univariate
logf <- function(x) log(x)
N <- 500
x <- exp(runif(N, 0, 2*pi))
y <- sin(log(x)) + rnorm(N, sd=.5)
xtrue <- exp(seq(0,2*pi,length=1000))
ytrue <- sin(log(xtrue))
mod.dt1 <- gam(y ~ s(x))
mod.dt2 <- gam(y ~ s(x, bs="dt", xt=list(tf=logf)))
mod.dt3 <- gam(y ~ s(x, bs="ad"))
fhat1 <- predict(mod.dt1, newdata=data.frame(x=xtrue))
fhat2 <- predict(mod.dt2, newdata=data.frame(x=xtrue))
fhat3 <- predict(mod.dt3, newdata=data.frame(x=xtrue))

plot(xtrue,ytrue,lwd=3, type="l", log="x", xlab="x", ylab="y")
lines(xtrue,fhat1,col="red")
lines(xtrue,fhat2,col="blue")
lines(xtrue,fhat3,col="green")
legend("topright", c("Truth", "f(x)", "f(log(x))", "f(x) adaptive"),
       col=c("black", "red", "blue", "green"), lty="solid", lwd=c(3,1,1,1))

N <- 500
x1 <- rnorm(N)
x2 <- rnorm(N, sd=2)
y <- sin(x1) + cos(x2) + rnorm(N, sd=.5)

sp1 <- s(x1)
sp2 <- s(x2)
mod.ms1 <- gam(y ~ s(x1, x2, bs="ms", xt="ps"))
mod.ms2 <- gam(y ~ s(x1) + s(x2))


y <- sin(2*pi*x/6) + rnorm(N, sd=.5)
mod <- gam(y ~ s(x, bs="nb"))

