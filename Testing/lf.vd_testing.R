library(devtools)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

dev_mode()
load_all("../refund")
load_all("../mgcvTrans/")

plotMe <- function(est, lims=NULL) {
  if (is.null(lims)) lims <- range(est$value)
  est$value[est$value<lims[1]] <- lims[1]
  est$value[est$value>lims[2]] <- lims[2]
  ggplot(est, aes(SOFA.arg, SOFA.vd)) +
    geom_tile(aes(colour=value, fill=value)) +
    scale_fill_gradientn(  name="", limits=lims, colours=rev(brewer.pal(11,"Spectral"))) +
    scale_colour_gradientn(name="", limits=lims, colours=rev(brewer.pal(11,"Spectral"))) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    #labs(title=as.character(substitute(est))) +
    theme_bw()
}



data(sofa)

# Untransformed
#TPRS
fit.vd1.1 <- pfr(death ~ lf.vd(SOFA) + age + los, family="binomial", data=sofa)
est.vd1.1 <- coef(fit.vd1.1, n=173, n2=173) %>% filter(SOFA.arg <= SOFA.vd)
plotMe(est.vd1.1, lims = c(-2,6))

# TPBS
fit.vd1.2 <- pfr(death ~ lf.vd(SOFA, basistype="te", bs="ps") + age + los, family="binomial", data=sofa)
est.vd1.2 <- coef(fit.vd1.2, n=173, n2=173) %>% filter(SOFA.arg <= SOFA.vd)
plotMe(est.vd1.2, lims = c(-2,6))


# Lagged
# TPRS
fit.vd2.1 <- pfr(death ~ lf.vd(SOFA, transform = "lagged") + age + los,
                 family="binomial", data=sofa)
est.vd2.1 <- coef(fit.vd2.1, n=173, n2=173) %>% filter(SOFA.arg <= SOFA.vd)
plotMe(est.vd2.1, c(-2,6))

#TPBS
fit.vd2.2 <- pfr(death ~ lf.vd(SOFA, transform = "lagged", basistype="te", bs="ps") + age + los,
                 family="binomial", data=sofa)
est.vd2.2 <- coef(fit.vd2.2, n=173, n2=173) %>% filter(SOFA.arg <= SOFA.vd)
plotMe(est.vd2.2, c(-2,6))


# Domain-Standardized
# TPRS:
fit.vd3.1 <- pfr(death ~ lf.vd(SOFA, transform = "standardized") + age + los,
                 family="binomial", data=sofa)
est.vd3.1 <- coef(fit.vd3.1, n=173, n2=173) %>% filter(SOFA.arg <= SOFA.vd)
plotMe(est.vd3.1, c(-2,4))

# TPBS: DOESN'T WORK
fit.vd3.2 <- pfr(death ~ lf.vd(SOFA, transform = "standardized", bs="ps", basistype = "te")
                 + age + los, family="binomial", data=sofa)
est.vd3.2 <- coef(fit.vd3.2, n=173, n2=173) %>% filter(SOFA.arg <= SOFA.vd)
plotMe(est.vd3.2, c(-2,4))



# No Interaction
fit.vd4.1 <- pfr(death ~ lf.vd(SOFA, transform = "noInteraction", bs="ps") + age + los,
                 family="binomial", data=sofa)
est.vd4.1 <- coef(fit.vd4.1, n=173, n2=173) %>% filter(SOFA.arg <= SOFA.vd)
plotMe(est.vd4.1)


# Linear Interaction
# One penalty
fit.vd5.1 <- pfr(death ~ lf.vd(SOFA, transform = "linear", bs="ps") + age + los,
                 family="binomial", data=sofa)
est.vd5.1 <- coef(fit.vd5.1, n=173, n2=173) %>% filter(SOFA.arg <= SOFA.vd)
plotMe(est.vd5.1, c(-2,6))

# Multiple penalties
fit.vd5.2 <- pfr(death ~ lf.vd(SOFA, transform = "linear", bs="ps", msp=TRUE)
                 + age + los, family="binomial", data=sofa)
est.vd5.2 <- coef(fit.vd5.2, n=173, n2=173) %>% filter(SOFA.arg <= SOFA.vd)
plotMe(est.vd5.2, c(-2,6))



# Quadratic Interaction
fit.vd6.1 <- pfr(death ~ lf.vd(SOFA, transform = "quadratic", bs="ps") + age + los,
                 family="binomial", data=sofa)
est.vd6.1 <- coef(fit.vd6.1, n=173, n2=173) %>% filter(SOFA.arg <= SOFA.vd)
plotMe(est.vd6.1, c(-2,6))

