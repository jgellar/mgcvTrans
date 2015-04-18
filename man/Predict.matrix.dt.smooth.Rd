\name{Predict.matrix.dt.smooth}
\alias{Predict.matrix.dt.smooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Predict.matrix.dt.smooth(object, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
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
function (object, data) 
{
    tf <- object$tf
    tdata <- lapply(tf, function(f) {
        args <- formals(f)
        argnms <- names(args)
        if (!(all(argnms) \%in\% object$term)) 
            argnms <- object$term[1:length(args)]
        if (!all(argnms \%in\% names(data))) {
            miss <- argnms[!(argnms \%in\% names(data))]
            stop(paste0("Variable(s) ", paste(argnms[miss], collapse = ", "), 
                " needed but not supplied to Predict.matrix"))
        }
        else {
            calldat <- data[argnms]
        }
        names(calldat) <- names(args)
        do.call(f, calldat)
    })
    untr <- names(data)[!(names(data) \%in\% names(tdata))]
    tdata[untr] <- data[untr]
    Predict.matrix(object, tdata)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
