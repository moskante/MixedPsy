\name{MERdelta.probit}
\alias{MERdelta.probit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MERdelta.probit(xplode.obj, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xplode.obj}{
%%     ~~Describe \code{xplode.obj} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
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
function (xplode.obj, alpha = 0.05) 
{
    if (xplode.obj$family$link != "probit") {
        output = NA
        print("Use a probit link function")
    }
    else {
        n.pf = length(xplode.obj$psychometrics)
        output = vector("list", length = n.pf)
        names(output) = names(xplode.obj$psychometrics)
        for (i in 1:n.pf) {
            pse <- -(xplode.obj$psychometrics[[i]]$intercept[1]/xplode.obj$psychometrics[[i]]$slope[1])
            slope <- xplode.obj$psychometrics[[i]]$slope[1]
            var.intercept <- xplode.obj$psychometrics[[i]]$intercept[2]
            var.slope <- xplode.obj$psychometrics[[i]]$slope[2]
            cov.intercept.slope <- xplode.obj$psychometrics$pf1$cov
            var.pse <- (1/slope^2) * (var.intercept + (2 * pse * 
                cov.intercept.slope) + (pse^2 * var.slope))
            inferior.pse <- pse - (qnorm(1 - (alpha/2)) * sqrt(var.pse))
            superior.pse <- pse + (qnorm(1 - (alpha/2)) * sqrt(var.pse))
            jnd <- qnorm(0.75) * (1/slope)
            var.jnd <- (qnorm(0.75) * (-1/slope^2))^2 * var.slope
            inferior.jnd <- jnd - (qnorm(1 - (alpha/2)) * sqrt(var.jnd))
            superior.jnd <- jnd + (qnorm(1 - (alpha/2)) * sqrt(var.jnd))
            output[[i]] <- matrix(rbind(c(pse, sqrt(var.pse), 
                inferior.pse, superior.pse), c(jnd, sqrt(var.jnd), 
                inferior.jnd, superior.jnd)), nrow = 2, dimnames = list(param <- c("pse", 
                "jnd"), statistics = c("Estimate", "Std. Error", 
                "Inferior", "Superior")))
        }
    }
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
