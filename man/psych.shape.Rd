\name{psych.shape}
\alias{psych.shape}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting and Plotting Psychometric Functions
}
\description{
\code{psych.shape()} is for plotting a psychometric function of a known pse and jnd on an existing plot.
}
\usage{
psych.shape(pse = 0, jnd, x.range = c(NA, NA),  ps.link = "probit", ps.col = "black", ps.lwd = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pse, jnd}{
the pse and the jnd of the desired psychometric function
}
  \item{ps.link}{
a link function for the binomial family of error distribution. See `Details'
}
  \item{x.range}{
a vector of length two specifying the range of the function
}
  \item{ps.col}{
color of the line to be plotted
}
  \item{ps.lty}{
line type the line to be plotted
}
}
\details{
psych.shape() can be used to visualize the predicted results of a psychophysical experiment or to plot a fitted psychometric function whose values of pse and jnd are known. Currently only working with probit and logit link function.
}
\value{
}
\references{
Moscatelli A, Mezzetti M, Lacquaniti F (2012). Modelling Psychophysical Data at the Population-Level: The Generalized Linear Mixed Model. Journal of Vision, 12(11):26, 1-17.
}
\author{
Alessandro Moscatelli, Priscilla Balestrucci
}
\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link[stats]{glm}} for for Generalized Linear Models.
}
\examples{
#Example
y = c(0,1)
x = c(-40, 40)
plot(y ~ x, type = "n", bty = "n", lab = c(5,3,7))
psych.shape(pse = 0, jnd = 6, x.range = c(-40, 40), ps.col = "gray", ps.lwd = 3)
psych.shape(pse = 6, jnd = 6, x.range = c(-40, 40), ps.col = "black")
psych.shape(pse = 6, jnd = 6, x.range = c(-40, 40), ps.col = "red", ps.link = "logit", ps.lwd = 3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
