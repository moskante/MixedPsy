\name{psych.function}
\alias{psych.function}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting and Plotting Psychometric Functions
}
\description{
\code{psych.function} is used to fit psychometric functions using either \code{glm()} or \code{brglm()}, estimate the PSE, the JND and the related confidence intervals and draw the curve on an existing plot.
}
\usage{
psych.function(ps.formula, ps.link, ps.data, x.range = c(NA, NA), ps.x = NA, ps.lines = F, ps.col = "black", ps.lty = "solid", br = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ps.formula}{
an object of class ``formula'', such as \code{cbind(yes, no) ~ X} 
}
  \item{ps.link}{
a link function for the binomial family of error distribution. See `Details'
}
  \item{ps.data}{
  a data frame including the variables in the model
}
  \item{ps.lines}{
logical. If TRUE, model predictions and confidence intervals of the PSE will be added to an existing plot.
}
\item{ps.x}{
optionally, a data frame in which to look for variables with which to predict. See `Details'
}
  \item{x.range}{
a vector of length two specifying the range for model predictions.
}
  \item{ps.col}{
color of the lines to be plotted
}
  \item{ps.lty}{
line type of the lines to be plotted
}
  \item{br}{
if TRUE, brglm is used if fitted values are equal to 0 or 1 
}

}
\details{
If lines is TRUE, the function draw model predictions on an existing plot. Only for univariable glm of the type \code{F(Y) ~ X}, where X is a continuous predictor. If ps.x is empty, the new data frame is a vector of length = 1000, whose range is specified from x.range. Std. Errors and 95\% confidence intervals of the PSE and JND are estimated via Delta Methods, see Faraggi et al (2003).
}
\value{
\code{psych.function} returns list including the fitted glm (or brglm), the estimate of PSE and JND and a flag to indicate if brglm was called.
}
\references{
Faraggi D, Izikson P, Reiser B (2003). Confidence intervals for the 50 per cent response dose. Statistics in Medicine, 22(12), 1977–1988.

Moscatelli A, Mezzetti M, Lacquaniti F (2012). Modelling Psychophysical Data at the Population-Level: The Generalized Linear Mixed Model. Journal of Vision, 12(11):26, 1-17.
}
\author{
Alessandro Moscatelli
}
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link[stats]{glm}} for for Generalized Linear Models.
}
\examples{
#this is to simulate data from a single participant
datafr.S1 <- MERsimulate(fixeff = c(-7.5, 0.0875), nsubject = 1, constant = T)

#fit a glm (probit link)
model.glm = glm(formula = cbind(Longer, Total - Longer) ~ X, family = binomial(link = "probit"), data = datafr.S1)

#psychometric function to fit single-subject data
plot(Longer/Total ~ X, data = datafr.S1)
fit.S1 = psych.function(ps.formula = cbind(Longer, Total - Longer) ~ X,
                        ps.link = "probit", ps.data = datafr.S1,
                        x.range = c(40, 120), ps.lines = T)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
