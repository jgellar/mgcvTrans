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


