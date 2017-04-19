\name{delta.psy.probit}
\alias{delta.psy.probit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PSE and JND  Using Delta Method
}
\description{
Estimates the Point of Subjective Equivalence (PSE), the Just Noticeable Difference (JND) and the related Standard Errors by means of Delta Method. The method apply to both, univariable GLM and GLMM having a \emph{probit} link function.
}
\usage{
delta.psy.probit(model, alpha = 0.05, lme4 = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
the fitted model. Either an object of class \code{"glm"} if \code{lme4 = F}  or an object of class \code{"\linkS4class{merMod}"} otherwhise. 
}
  \item{alpha}{
significance level of the confidence interval.
}
  \item{lme4}{
a logical indicating whether fitted model is an object of class \code{"\linkS4class{merMod}"}
}
}
\details{
\code{delta.psy.probit} estimates PSE and JND of a univariable psychometric function (object of class \code{"glm"}) or of a univarible GLMM (object of class \code{"\linkS4class{merMod}"}). Use \code{MERtreatment} for multivariable GLMM.
}
\value{
%%  ~Describe the value returned
\code{delta.psy.probit} returns a matrix including Estimate, Standard Error, Inferior and Superior Confidence Interval of PSE and JND. Confidence Intervals are computed as: \eqn{Estimate +/- z(1-(\alpha/2)) * Std.Error}.
}
\references{
Faraggi D, Izikson P, Reiser B (2003). Confidence intervals for the 50 per cent response dose. Statistics in Medicine, 22(12), 1977–1988.

Moscatelli A, Mezzetti M, Lacquaniti F (2012). Modelling Psychophysical Data at the Population-Level: The Generalized Linear Mixed Model. Journal of Vision, 12(11):26, 1-17.
}
\author{
Alessandro Moscatelli
}
\note{
The function assumes that the first model coefficient is the intercept and the second is the slope. The estimate of the JND assumes a \emph{probit} link function.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link[stats]{glm}} for for Generalized Linear Models (without random effects) and \code{\link[lme4]{glmr}} for Generalized Linear Mixed Models (including random effects). \code{\link{pseMer}} provides the bootstrap-based confidence intervals.
}
\examples{
#this is to simulate data from a single participant
datafr.S1 <- MERsimulate(fixeff = c(-7.5, 0.0875), nsubject = 1, constant = T)
#fit a glm (probit link)
model.glm = glm(formula = cbind(Longer, Total - Longer) ~ X, family = binomial(link = "probit"), data = datafr.S1)
delta.psy.probit(model.glm)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
