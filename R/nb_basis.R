#' Smooths with No Basis
#' 
#' The \code{nb} "basis" uses a separate coefficient for each value of
#' the variables this is a smooth function of. This is equivalent to using
#' an identity matrix as the basis matrix. The coefficients can be penalized
#' with a difference penalty. This basis is most useful when there are not very
#' many unique values of the variable, making dimension reduction unnecessary.
#' 
#' 
#' 
#' 
#' 
#' @examples
#' N <- 100
#' x <- sample(-3:3, N, replace=T)
#' y <- sin(2*pi*x/6) + rnorm(N, sd=.5)
#' mod <- gam(y ~ s(x, bs="nb"))
#' 

smooth.construct.nb.smooth.spec <- function(object, data, knots) {
  x <- factor(data[[1]])
  m <- object$p.order
  k <- object$bs.dim
  if (k!=-1 & k!=length(levels(x)))
    warning(paste0("setting k to the number of unqiue values of the variable",
                   "to be smoothed."))
  k <- length(levels(x))
  
  # Difference Penalty
  if (is.na(m)) m <- 2 # Default
  L <- diag(k)
  for (i in 1:m) L <- diff(L)
  L1 <- L[nrow(L),]
  for (i in 1:m) {
    L1 <- c(0, L1[-k])
    L <- rbind(L, L1)
  }
  D <- t(L) %*% L
  X <- model.matrix(~x)
  #colnames(X) <- paste0("s(", object$term[1], ").", 1:k)
  
  # Modify return object
  object$X <- X
  object$S <- list(D)
  object$rank <- qr(D)$rank
  object$null.space.dim <- k - object$rank
  object$df <- k #need this for gamm
  object$argvals <- as.numeric(levels(x))
  object$interpolation <- "linear"
  class(object) <- "nb.smooth"
  object
}

#' Predict.matrix method for nb basis

Predict.matrix.nb.smooth <- function(object, data) {
  # Prediction method for parameteric nb basis
  if (all(data[[1]]) %in% object$argvals) {
    # No interpolation needed
    newx <- factor(data[[1]], levels=object$argvals)
    model.matrix(~newx)
  } else {
    # Requires interpolation!
    f0 <- model.matrix(~factor(object$argvals))
    sapply(1:nrow(f0), function(j) {
      if (object$interpolation=="linear") {
        approx(object$argvals,f0[,j], data[[1]])$y
      } else if (object$interpolation %in% c("loess", "lowess")) {
        predict(loess(y ~ x, data=data.frame(y=f0[,j], x=object$argvals)),
                newdata=data.frame(x=data[[1]]))
      } else {
        stop(paste0("Interpolation method ", object$interpolation,
                    " not supported."))
      }
    })
#     
#     if (object$interpolation=="linear") {
#       W <- t(sapply(sapply(data[[1]], function(x) {
#         approx(object$argvals, 1:nrow(f0), x)$y
#       }), function(y) {
#         d <- abs(sapply(1:nrow(f0), function(z) {
#           y-z
#         }))
#         d[d>1] <- 1
#         1-d
#       }))
#       W %*% f0
#     } else if (object$interpolation=="lowess") {
#       sapply(1:nrow(f0), function(j) {
#         lfit <- loess(y ~ x, data=data.frame(y=f0[,j], x=object$argvals))
#         tmp=predict(lfit, newdata=data.frame(x=data[[1]]))
#       })
#     } else {
#       stop("Interpolation method not supported")
#     }
  }
}