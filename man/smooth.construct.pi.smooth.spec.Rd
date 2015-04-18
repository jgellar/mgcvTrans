\name{smooth.construct.pi.smooth.spec}
\alias{smooth.construct.pi.smooth.spec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
smooth.construct.pi.smooth.spec(object, data, knots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{knots}{
%%     ~~Describe \code{knots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, data, knots) 
{
    if (length(object$term) < 2) 
        stop("\"pi\" smooth must have at least two variables")
    xt <- object$xt
    bs = (newxt <- NULL)
    idx <- length(object$term)
    if (is.character(xt)) {
        tf <- xt
    }
    else if (is.list(xt)) {
        if (is.null(xt$tf)) {
            tf <- xt
        }
        else {
            if (!all(names(xt) \%in\% c("tf", "bs", "xt", "idx"))) 
                warning("Ignoring unrecognized elements of \"xt\"")
            tf <- xt$tf
            bs <- xt$bs
            newxt <- xt$xt
            if (!is.null(xt$idx)) 
                idx <- xt$idx
        }
    }
    else stop("xt must be a character string or list")
    if (is.character(tf)) {
        if (tf == "none") 
            tf <- list(function(t) t^0)
        else if (tf == "linear") 
            tf <- list(function(t) t^0, function(t) t^1)
        else if (tf == "quadratic") 
            tf <- list(function(t) t^0, function(t) t^1, function(t) t^2)
        else stop("Unrecognized character string for transformation functions")
    }
    else if (is.list(tf)) {
        if (any(sapply(tf, class) != "function")) 
            stop("Incorrect input for transformation functions")
    }
    else stop("Incorrect input for transformation functions")
    tvar <- data[[object$term[idx]]]
    dat <- data[object$term[-idx]]
    n_transform <- length(tf)
    tf_name <- paste0(object$term[idx], "_tf_", 1:n_transform)
    t_X <- vapply(tf, function(f) do.call(f, list(t = tvar)), 
        numeric(length(tvar)))
    args <- list(k = object$bs.dim, fx = object$fixed, m = object$p.order, 
        xt = newxt, id = object$id, sp = object$sp)
    if (!is.null(bs)) 
        args$bs <- bs
    if (object$by != "NA") {
        args$by <- data[[object$by]]
        dat[[object$by]] <- data[[object$by]]
    }
    callargs <- sapply(object$term[-idx], as.name)
    names(callargs) <- NULL
    smoothspec <- do.call(mgcv::s, append(callargs, args))
    sm0 <- smooth.construct(smoothspec, data = dat, knots = knots)
    sm <- sm0
    sm$term <- c(object$term[-idx], object$term[idx])
    sm$bs.dim <- sm$bs.dim * n_transform
    sm$null.space.dim <- sm$null.space.dim * n_transform
    sm$df <- sm$df * n_transform
    sm$dim <- length(object$term)
    sm$label <- paste0("f(", object$term[idx], ")*", sm$label)
    sm$xt <- object$xt
    sm$X <- mgcv::tensor.prod.model.matrix(list(t_X, sm$X))
    sm$S[[1]] <- diag(n_transform) \%x\% sm$S[[1]]
    sm$bs <- strsplit(class(smoothspec), ".", fixed = T)[[1]][1]
    sm$tf <- tf
    sm$sm <- sm0
    if (!is.null(sm$drop.null)) {
        if (sm$drop.null > 0) {
            ind <- 1:(sm$bs.dim - sm$null.space.dim)
            if (FALSE) {
                np <- nat.param(sm$X, sm$S[[1]], rank = sm$bs.dim - 
                  sm$null.space.dim, type = 0)
                sm$P <- np$P
                sm$S[[1]] <- diag(np$D)
                sm$X <- np$X[, ind]
            }
            else {
                sm$S[[1]] <- sm$S[[1]][ind, ind]
                sm$X <- sm$X[, ind]
                sm$cmX <- colMeans(object$X)
                sm$X <- sweep(object$X, 2, object$cmX)
            }
            sm$null.space.dim <- 0
            sm$df <- sm$df - M
            sm$bs.dim <- sm$bs.dim - M
            sm$C <- matrix(0, 0, ncol(object$X))
        }
    }
    class(sm) <- "pi.smooth"
    sm
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
