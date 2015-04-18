\name{smooth.construct.dt.smooth.spec}
\alias{smooth.construct.dt.smooth.spec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
smooth.construct.dt.smooth.spec(object, data, knots)
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
    tf <- if (class(object) \%in\% c("tensor.smooth.spec", "t2.smooth.spec")) {
        x <- sapply(chk$margin, function(mar) {
            mar$xt$tf
        })
        x[!sapply(x, is.null)][[1]]
    }
    else {
        object$xt
    }
    nt <- length(object$term)
    if (!is.list(tf)) 
        tf <- list(tf)
    tf <- lapply(tf, function(f) {
        if (is.character(f) & nt == 2) {
            if (f == "s-t") {
                function(s, t) s - t
            }
            else if (f == "s/t") {
                function(s, t) s/t
            }
            else {
                stop("Unrecognized character string for 'dt' transformation function")
            }
        }
        else if (is.function(f)) {
            f
        }
        else stop("Unrecognized type for 'dt' transformation function")
    })
    if (is.null(names(tf))) 
        names(tf) <- object$term[1:length(tf)]
    tdata <- lapply(tf, function(f) {
        args <- formals(f)
        argnms <- names(args)
        if (length(args) > length(object$term)) 
            stop("Transformation function of too many arguments is supplied")
        if (!(all(argnms) \%in\% object$term)) 
            argnms <- object$term[1:length(args)]
        if (!all(argnms \%in\% names(data))) {
            miss <- argnms[!(argnms \%in\% names(data))]
            stop(paste0("Variable(s) ", paste(argnms[miss], collapse = ", "), 
                " needed but not supplied to smooth constructor"))
        }
        else {
            calldat <- data[argnms]
        }
        names(calldat) <- names(args)
        do.call(f, calldat)
    })
    untr <- names(data)[!(names(data) \%in\% names(tdata))]
    tdata[untr] <- data[untr]
    if (class(object) \%in\% c("tensor.smooth.spec", "t2.smooth.spec")) {
        object$margin <- lapply(object$margin, function(mar) {
            xt <- mar$xt
            bs <- ifelse(is.null(xt$bs), "tp", xt$bs)
            class(mar) <- paste0(bs, ".smooth.spec")
            mar$xt <- xt[!(names(xt) \%in\% c("tf", "bs"))]
            mar
        })
    }
    else {
        xt <- object$xt
        bs <- ifelse(is.null(xt$bs), "tp", xt$bs)
        class(object) <- paste0(bs, ".smooth.spec")
        object$xt <- xt[!(names(xt) \%in\% c("tf", "bs"))]
    }
    sm <- smooth.construct(object, data = tdata, knots = knots)
    sm$tf <- tf
    sm$class <- class(sm)
    class(sm) <- "dt.smooth"
    sm
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
