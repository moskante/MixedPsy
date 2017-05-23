% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/psych_utils.R
\name{PsychShape}
\alias{PsychShape}
\title{Fitting and Plotting Psychometric Functions}
\usage{
PsychShape(pse = 0, jnd, x.range = c(NA, NA), ps.link = "probit",
  ps.col = "black", ps.lwd = 1)
}
\arguments{
\item{pse, jnd}{the pse and the jnd of the desired psychometric function}

\item{x.range}{a vector of length two specifying the range of the function}

\item{ps.link}{a link function for the binomial family of error distribution.
See `Details'}

\item{ps.col}{color of the line to be plotted}

\item{ps.lty}{line type the line to be plotted}
}
\description{
\code{PsychShape()} plots a psychometric function of a known pse and jnd
on an existing plot.
}
\details{
PsychShape() can be used to visualize the predicted results of a
psychophysical experiment or to plot a fitted psychometric function whose
values of pse and jnd are known. Currently only working with probit and logit
link function.
}
\examples{
y = c(0,1)
x = c(-40, 40)
plot(y ~ x, type = "n", bty = "n", lab = c(5,3,7))
PsychShape(pse = 0, jnd = 6, x.range = c(-40, 40), ps.col = "gray", ps.lwd = 3)
PsychShape(pse = 6, jnd = 6, x.range = c(-40, 40), ps.col = "black")
PsychShape(pse = 6, jnd = 6, x.range = c(-40, 40), ps.col = "red", ps.link = "logit", ps.lwd = 3)

}
\references{
Moscatelli A, Mezzetti M, Lacquaniti F (2012). Modelling Psychophysical Data
at the Population-Level: The Generalized Linear Mixed Model.
Journal of Vision, 12(11):26, 1-17.
}
\seealso{
\code{\link[stats]{glm}} for for Generalized Linear Models.
}