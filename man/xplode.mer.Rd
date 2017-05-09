\name{xplode.mer}
\alias{xplode.mer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract values from a fitted GLMM objects
}
\description{
A function to extract values from an object of class \code{"\linkS4class{merMod}"} (more specifically, from an object of subclass glmerMod).
}
\usage{
xplode.mer(model, name.cont = NA, name.factor = NA, names.response = NA,
            define.pf = list(pf1 = list(intercept = 1, slope = 2)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
the GLMM fitted with \code{glmer}. An object of class \code{"\linkS4class{merMod}"}.
}
  \item{name.cont}{
a string providing the name of the continuous predictor, as in the formula object of the fitted model
}
  \item{name.factor}{
  a string providing the name of name of the categorical predictor, as in the formula object of the fitted model
}
  \item{names.response}{
optional. A string providing the name of name of the response variable, as in the formula object of the fitted model
}
  \item{define.pf}{
optional. Specifies which parameter pertains to the intercept and which to the slope in the formula object.
}
}
\details{
For simplicity, several \code{MixedPsy} functions take as input an object of class \code{xplode} instead of an object of class \code{"\linkS4class{merMod}"}. Most of these functions assume by default that the continuous predictor is entered first in the formula object. It is possible to use a different order, this requires to specify which parameter pertains to the intercept and which to the slope, by changing \code{define.pf}.
}
\value{
An object of class \code{"\linkS4class{merMod}"} to be used with other \code{MixedPsy} functions.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Alessandro Moscatelli, Priscilla Balestrucci
}
\note{
%%  ~~further notes~~
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
datafr = MERsimulate(nsubjects = 10)
mod1 = glmer(formula = cbind(Longer, Total - Longer) ~ X + (1 | Subject), family = binomial(link = "probit"), data = datafr)
xplode.mod1 = xplode.mer(model = mod1, name.cont = "X")
MERdelta.probit(xplode.mod2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
