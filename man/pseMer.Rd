\name{pseMer}
\alias{pseMer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pseMer(mer.obj, B = 200, FUN = NULL, beep = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mer.obj}{
%%     ~~Describe \code{mer.obj} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{FUN}{
%%     ~~Describe \code{FUN} here~~
}
  \item{beep}{
%%     ~~Describe \code{beep} here~~
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
function (mer.obj, B = 200, FUN = NULL, beep = F)
{
    boot.flag = require("boot")
    if (boot.flag == F) {
        print("Package boot is missing. Install the package to use the function.")
    }
    if (is.null(FUN)) {
        myfun <- function(mer.obj) {
            jndpse = vector(mode = "numeric", length = 2)
            names(jndpse) = c("JND", "PSE")
            jndpse[1] = qnorm(0.75)/fixef(mer.obj)[2]
            jndpse[2] = -fixef(mer.obj)[1]/fixef(mer.obj)[2]
            return(jndpse)
        }
    }
    else {
        myfun = match.fun(FUN)
    }
    np = length(myfun(mer.obj))
    parname = names(myfun(mer.obj))
    if (is.null(parname)) {
        print("Warning messages: Parameters have no names")
    }
    summary = matrix(NA, nrow = np, ncol = 3, dimnames = list(parname,
        c("Estimate", "Inferior", "Superior")))
    boot.samp <- bootMer(mer.obj, myfun, nsim = B)
    summary[, 1] = boot.samp$t0
    jndpseconf = vector(mode = "list", length = np)
    for (i in 1:np) {
        jndpseconf[[i]] <- boot.ci(boot.samp, type = c("norm",
            "basic", "perc"), index = i)
        print(paste(parname[i], " 95\% CI:", jndpseconf[[i]]$perc[4],
            "  ", jndpseconf[[i]]$perc[5]))
        summary[i, 2] = jndpseconf[[i]]$perc[4]
        summary[i, 3] = jndpseconf[[i]]$perc[5]
    }
    if (beep == T) {
        beep.flag = require("beepr")
        if (beep.flag == F) {
            print("Package beep is missing. Install the package or set beep = F")
        }
        beep()
    }
    out = list(summary, boot.samp, jndpseconf)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
