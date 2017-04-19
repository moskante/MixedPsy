\name{kombo}
\alias{kombo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
kombo(vector)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vector}{
%%     ~~Describe \code{vector} here~~
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
function (vector) 
{
    n = length(vector)
    if (n > 1) {
        N = length(2:n)
        temp = vector("list", N)
        twice = matrix(NA, nrow = 2, ncol = n)
        for (i in 1:n) {
            twice[, i] = rep(vector[i], 2)
        }
        temp[[1]] = twice
        for (j in 2:n) {
            temp[[j]] = combn(vector, j)
        }
    }
    else {
        temp = vector
    }
    return(temp)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
