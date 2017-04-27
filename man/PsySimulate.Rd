\name{Psysimulate}
\alias{PsySimulate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate psychophysical data.
}
\description{
Given the arrays of fixed and random effects, as well as the covariance, and the characteristic of the simulated experiment (i.e., ) the function simulates a dataset in which for each subject the following information is provided: the slope and intercept value of the psychometric function, and the simulated responses to the stimulus levels that fit that function.
}
\usage{
MERsimulate(fixeff, raneff, nsubjects, pps, ntrials, xint, constant = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fixeff}{
Array of fixed effects. First item is the intercept, second item is the slope. 
}
  \item{raneff}{
Array of random effects. First item is the intercept, second item is the covariance, third item is the slope. 
}
  \item{nsubjects}{
Number of subjects to simulate data for. Default is 8. 
}
  \item{ntrials}{
Number of trials for each stimulus level. Default is 40.
}
  \item{xint}{
Range of the stimulus interval. Default is c(40,120).
}
  \item{pps}{
Number of stimulus levels. Default is 9
}
  \item{constant}{If set to FALSE, stimulus levels are randomly generated, uniformly distributed values within the selected interval.
  If constant = TRUE, the X interval is divided in  intervals of constant length. Default is FALSE.
%%     ~~Describe \code{constant} here~~
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


#simulate dataset (one subject)
datafr.S1 <- MERsimulate(nsubject = 1, constant = T)     
 
#psychometric function to fit single-subject data
plot(Longer/Total ~ X, data = datafr.S1)
fit.S1 = psych.function(ps.formula = cbind(Longer, Total - Longer) ~ X,
                        ps.link = "probit", ps.data = datafr.S1,
                        x.range = c(40, 120), ps.lines = T)
fit.S1 


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
