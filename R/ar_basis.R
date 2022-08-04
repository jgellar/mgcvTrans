#' Autoregressive Penalty
#' 
#' @export
#' @importFrom stats approx
#' @importFrom stats loess
#' 

smooth.construct.ar.smooth.spec <- function(object, data, knots) {
  k <- object$bs.dim
  interp <- ifelse(is.null(object$xt), "linear", object$xt)
  knots <- knots[[object$term]]
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
        predict(mgcv::gam(y ~ s(x, bs="ps"), data=data.frame(y=f0[,j], x=knots)),
                newdata=data.frame(x=data[[1]]))
        #stop("b-spline interpolation not yet implemented")
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
  class(object) <- "ar.smooth"
  object
}

#' Predict.matrix method for ar basis
#' 
#' @param object a \code{ar.smooth} object created by
#'   \code{\link{smooth.construct.ar.smooth.spec}}, see
#'   \code{\link[mgcv]{smooth.construct}}
#' @param data  see \code{\link[mgcv]{smooth.construct}}
#' @return design matrix for \code{dt} terms
#' @author Jonathan Gellar
#' @export
#' @importFrom stats model.matrix approx predict
#' 
Predict.matrix.ar.smooth <- function(object, data) {
  # Prediction method for parameteric ar basis
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
      } else if (object$interpolation %in% c("bs", "bspline")) {
        predict(mgcv::gam(y ~ s(x, bs="ps", k=object$bs.dim, sp=0),
                    data=data.frame(y=f0[,j], x=object$knots)),
                newdata=data.frame(x=data[[1]]))
      } else {
        stop(paste0("Interpolation method ", object$interpolation,
                    " not supported."))
      }
    })
  }
}
