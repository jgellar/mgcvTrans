smooth.construct.ms.smooth.spec <- function(object, data, knots) {
  obj1 = obj2 <- object
  xt <- object$xt
  
  smooths <- lapply(1:length(object$term), function(i) {
    obj.i <- object
    obj.i$term <- object$term[i]
    obj.i$dim <- 1
    obj.i$label <- paste0("s(", obj.i$term, ")")
    obj.i$xt <- NULL
    class(obj.i) <- paste0(xt, ".smooth.spec")
    dat.i <- data[obj.i$term]
    if (obj.i$by != "NA")
      dat.i[obj.i$by] <- data[obj.i$by]
    smooth.construct(obj.i, data=dat.i, knots=knots)
  })
  
  # Create new smooth object
  object$mar <- smooths
  object$X <- do.call("cbind", lapply(smooths, function(x) x$X))
  object$S <- lapply(1:length(smooths), function(k) {
    diag(as.numeric((1:length(smooths)) == k)) %x% smooths[[k]]$S[[1]]
  })
  object$bs.dim <- sum(sapply(smooths, function(x) x$bs.dim))
  object$rank <- sapply(object$S, function(D) qr(D)$rank)
  object$null.space.dim <- object$bs.dim - sum(object$rank)
  object$p.order <- sapply(smooths, function(x) x$p.order)
  class(object) <- "ms.smooth"
  object
}

