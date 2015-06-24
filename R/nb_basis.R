#' Smooths with No Basis
#' 
#' The \code{nb} "basis" uses a separate coefficient for each value of
#' the variables this is a smooth function of. The coefficients are penalized
#' with a difference penalty. This basis is most useful when there are not very
#' many unique values of the variable, making dimension reduction unnecessary.
#' 
#' @param object the object
#' @param data the data
#' @param knots the knots
#' 
#' 
#' @examples
#' 
#' \dontrun{
#' N <- 100
#' x <- sample(-3:3, N, replace=T)
#' y <- sin(2*pi*x/6) + rnorm(N, sd=.5)
#' mod <- gam(y ~ s(x, bs="nb"))
#' pre <- predict(mod, newdata=data.frame(x=c(-3:3)))
#' pre <- predict(mod, newdata=data.frame(x=seq(-3,3,length=100)))
#' mod <- gam(y ~ s(x, bs="nb", xt="lowess"))
#' pre <- predict(mod, newdata=data.frame(x=seq(-3,3,length=100)))
#' }
#' 
#' 

smooth.construct.nb.smooth.spec <- function(object, data, knots) {
  k <- object$bs.dim
  interp <- ifelse(is.null(object$xt), "linear", object$xt)
  if (is.null(knots)) {
    knots <- unique(data[[1]])
  }
  knots <- knots[order(knots)]
  
  if (k!=-1 & k!=length(knots)) {
    warning("ignoring k due to mismatch with the number of knots")
  }
  k <- length(knots)
  
  # Model matrix
  X <- if (all(data[[1]] %in% knots)) {
    # No interpolation needed
    model.matrix(~factor(data[[1]], levels=knots))
  } else {
    # Requires interpolation
    f0 <- model.matrix(~factor(knots))
    sapply(1:nrow(f0), function(j) {
      if (tolower(interp)=="linear") {
        approx(knots, f0[,j], data[[1]])$y
      } else if (tolower(interp) %in% c("loess", "lowess")) {
        predict(loess(y ~ x, data=data.frame(y=f0[,j], x=knots)),
                newdata=data.frame(x=data[[1]]))
      } else if (tolower(interp)%in% c("bs", "bspline")) {
        stop("b-spline interpolation not yet implemented")
      } else {
        stop(paste0("Interpolation method ", object$interpolation,
                    " not supported."))
      }
    })
  }
  
  # Difference Penalty
  m <- object$p.order
  if (is.na(m)) m <- 2 # Default
  L <- diag(k)
  for (i in 1:m) L <- diff(L)
  L1 <- L[nrow(L),]
  for (i in 1:m) {
    L1 <- c(0, L1[-k])
    L <- rbind(L, L1)
  }
  D <- t(L) %*% L
  
  # Modify return object
  object$X <- X
  object$S <- list(D)
  object$knots <- knots
  object$bs.dim <- k
  object$rank <- qr(D)$rank
  object$null.space.dim <- k - object$rank
  object$df <- k #need this for gamm
  object$interpolation <- interp
  class(object) <- "nb.smooth"
  object
}

#' Predict.matrix method for nb basis
#' 
#' @param object the object
#' @param data the data
#' 
Predict.matrix.nb.smooth <- function(object, data) {
  # Prediction method for parameteric nb basis
  if (all(data[[1]] %in% object$knots)) {
    # No interpolation needed
    newx <- factor(data[[1]], levels=object$knots)
    model.matrix(~newx)
  } else {
    # Requires interpolation!
    f0 <- model.matrix(~factor(object$knots))
    sapply(1:nrow(f0), function(j) {
      if (object$interpolation=="linear") {
        approx(object$knots,f0[,j], data[[1]])$y
      } else if (object$interpolation %in% c("loess", "lowess")) {
        predict(loess(y ~ x, data=data.frame(y=f0[,j], x=object$knots)),
                newdata=data.frame(x=data[[1]]))
      } else {
        stop(paste0("Interpolation method ", object$interpolation,
                    " not supported."))
      }
    })
  }
}



#     
#     if (object$interpolation=="linear") {
#       W <- t(sapply(sapply(data[[1]], function(x) {
#         approx(object$knots, 1:nrow(f0), x)$y
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
#         lfit <- loess(y ~ x, data=data.frame(y=f0[,j], x=object$knots))
#         tmp=predict(lfit, newdata=data.frame(x=data[[1]]))
#       })
#     } else {
#       stop("Interpolation method not supported")
#     }
