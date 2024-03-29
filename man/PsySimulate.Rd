% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PsySimulate.R
\name{PsySimulate}
\alias{PsySimulate}
\title{Simulate psychophysical data}
\usage{
PsySimulate(
  fixeff = c(-7, 0.0875),
  raneff = c(2.4, -0.002, 2e-06),
  nsubjects = 8,
  ntrials = 40,
  nintervals = 9,
  xint = c(40, 120),
  constant = T
)
}
\arguments{
\item{fixeff}{array of fixed effects. First item is the intercept, second
item is the slope.}

\item{raneff}{array of random effects. First item is the intercept, second
item is the covariance, third item is the slope.}

\item{nsubjects}{number of subjects to simulate. Default is 8.}

\item{ntrials}{number of trials for each stimulus level. Default is 40.}

\item{nintervals}{number of stimulus levels. Default is 9.}

\item{xint}{range of the stimulus interval. Default is c(40,120)}

\item{constant}{logical. If set to FALSE, stimulus levels are randomly generated,
uniformly distributed values within the selected interval.
otherwise, the X interval is divided in  intervals of constant
length. Default is TRUE.}
}
\value{
\code{PsySimulate} returns a simulated dataset. If no input arguments are specified, the function returns 
a dataset with the same characteristics as \code{\link{simul_data}}.
}
\description{
The function simulates data of a typical psychophysics experiment. For each simulated 
participant, the function returns the following information: individual slope and intercept 
coefficients, given the fixed and random effects parameters provided as input; summary of 
the simulated binomial response to a range of intensity levels between a specified range.
}
\examples{
datafr.S1 <- PsySimulate(fixeff = c(0, 1), xint = c(-5,5), 
nsubject = 1, ntrials = 60, nintervals = 10, constant = FALSE)
library(ggplot2)
g <- ggplot(datafr.S1, aes(X,Longer/Total)) + geom_point()

PsychShape(pse = 0, jnd = qnorm(0.75)/1, ps.link = "probit", 
x.range = c(-5,5), addTo = g, ps.color = "red")

}
\seealso{
\code{\link{PsychShape}} for plotting a psychometric function given PSE and JND.
}
