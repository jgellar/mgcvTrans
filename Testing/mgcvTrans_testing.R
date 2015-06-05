# Testing....

library(devtools)
dev_mode()
load_all()
#load_all("../pcox/")
load_all("../refund/")
#load_all("../FDA - Survival/pcox/")

library(mgcv)
library(survival)

# data(sofa)
# tmp <- lf(sofa$SOFA)
# cl <- tmp$call
# cl$bs <- "dt"
# 
# #cl[[5]] <- "dt"
# #names(cl)[4] <- "bs"
# 
# spec <- eval(cl, envir = tmp$data)
# sm <- smoothCon(spec, data = tmp$data, knots = NULL, absorb.cons = T)

data(DTI)
DTI1 <- DTI[DTI$visit==1 & complete.cases(DTI),]
f <- function(x) x^2

tmp <- lf(DTI1$cca, k=30, bs="dt", xt=list(tf=f))

fit.lf <- pfr(pasat ~ lf(cca, k=30, bs="dt", xt=list(tf=f)), data=DTI1)

chk <- te(var1, var2, xt=list(list(tf="s-t", bs="tp"), list(bs="cr")))
chk <- te(var1, var2, xt=list(bs=c("tp","cr")))


# dt testing - Univariate





y <- sin(2*pi*x/6) + rnorm(N, sd=.5)
mod <- gam(y ~ s(x, bs="nb"))






# Linear functional terms via gm




# Variable-domain funtional terms
data(sofa)
N <- nrow(sofa)
J <- ncol(sofa$SOFA)
J.i <- apply(sofa$SOFA, 1, function(x) {
  max(which(!is.na(x)))
})
tmat <- matrix(0:(J-1), nrow=N, ncol=J, byrow=T)/(J-1)
Tmat <- matrix((J.i-1), nrow=N, ncol=J)/(J-1)
Xmat <- sofa$SOFA
L  <- refundDevel:::getL(tmat, "simpson", J.i)
LX <- Xmat * L

tmat[is.na(Xmat)] <- tmat[!is.na(Xmat)][1]
Tmat[is.na(Xmat)] <- Tmat[!is.na(Xmat)][1]
LX[is.na(Xmat)] <- 0


fit.vd <- gam(death ~ s(tmat, Tmat, by=LX) + age + los,
              family="binomial", data=sofa)
est.vd <- coef.pfr(fit.vd, n=173, n2=173) %>% filter(tmat <= Tmat)

lims <- range(est.vd$value, na.rm=T)
lims <- c(-2,8) * 173
ggplot(est.vd, aes(tmat, Tmat)) +
  geom_tile(aes(fill=value, colour=value)) +
  scale_fill_gradientn(name="", limits=lims,
                       colours=rev(brewer.pal(11,"Spectral"))) +
  scale_colour_gradientn(name="", limits=lims,
                         colours=rev(brewer.pal(11,"Spectral"))) +
  theme_bw() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0))


fit.vd1 <- pfr(death ~ lf.vd(SOFA) + age + los, family="binomial",data=sofa)


lims <- c(-2,8)
ggplot(est.vd1, aes(SOFA.tmat, SOFA.Tmat)) +
  geom_tile(aes(fill=value, colour=value)) +
  scale_fill_gradientn(name="", limits=lims,
                       colours=rev(brewer.pal(11,"Spectral"))) +
  scale_colour_gradientn(name="", limits=lims,
                         colours=rev(brewer.pal(11,"Spectral"))) +
  theme_bw() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0))
