\name{psych.delta}
\alias{psych.delta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PSE and JND  Using Delta Method
}
\description{
Estimates the Point of Subjective Equivalence (PSE), the Just Noticeable Difference (JND) and the related Standard Errors by means of Delta Method. The method only applies to univariable GLMs (psychometric functions) having a \emph{probit} link function.
}
\usage{
psych.delta(model, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
the fitted psychometric function. An object of class \code{"glm"}. 
}
  \item{alpha}{
significance level of the confidence interval.
}
}
\details{
\code{psych.delta} estimates PSE and JND of a univariable psychometric function (object of class \code{"glm"}).
}
\value{
%%  ~Describe the value returned
\code{psych.delta} returns a matrix including Estimate, Standard Error, Inferior and Superior Confidence Interval of PSE and JND. Confidence Intervals are computed as: \eqn{Estimate +/- z(1-(\alpha/2)) * Std.Error}.
}
\references{
Faraggi D, Izikson P, Reiser B (2003). Confidence intervals for the 50 per cent response dose. Statistics in Medicine, 22(12), 1977–1988.

Moscatelli A, Mezzetti M, Lacquaniti F (2012). Modelling Psychophysical Data at the Population-Level: The Generalized Linear Mixed Model. Journal of Vision, 12(11):26, 1-17.
}
\author{
Alessandro Moscatelli, Priscilla Balestrucci
}
\note{
The function assumes that the first model coefficient is the intercept and the second is the slope. The estimate of the JND assumes a \emph{probit} link function.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link[stats]{glm}} for for Generalized Linear Models (without random effects) and \code{\link[lme4]{glmr}} for Generalized Linear Mixed Models (including random effects).
\code{MERdelta.probit} and \code{MERtreatment} for univarible and multivariable GLMM, respectively (object of class \code{"\linkS4class{merMod}"}).
\code{\link{pseMer}} provides the bootstrap-based confidence intervals.
}
\examples{
#this is to simulate data from a single participant
datafr.S1 <- MERsimulate(fixeff = c(-7.5, 0.0875), nsubject = 1, constant = T)
#fit a glm (probit link)
model.glm = glm(formula = cbind(Longer, Total - Longer) ~ X, family = binomial(link = "probit"), data = datafr.S1)
psych.delta(model.glm)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
