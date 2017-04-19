\name{MERtreatment}
\alias{MERtreatment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PSE/JND for Multivariable GLMM Using Delta Methods
}
\description{
Estimates the Point of Subjective Equivalence (PSE), the Just Noticeable Difference (JND) and the related Standard Errors by means of Delta Method. The method apply to multivariable GLMM having a \emph{probit} link function. The function is based on a recursive use of \code{glmer} and \code{delta.psy.probit}
}
\usage{
MERtreatment(xplode.obj, datafr, tr.alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xplode.obj}{
an object of class \code{xplode.obj}. The fitted model (object of class \code{"\linkS4class{merMod}"}) from \code{xplode.obj} includes one continuous predictor and one factorial predictor.
}
  \item{datafr}{
the data frame fitted with the GLMM model
}
  \item{tr.alpha}{
significance level of the confidence interval.
}
}
\details{
The function \code{MERtreatment} is based on a recursive use of \code{glmer} and \code{delta.psy.probit} to multivariable GLMM including continuous and factorial predictors. The same caveats of \code{delta.psy.probit} appliy (e.g., confidence interval based on normality assumption).
}
\value{
A list, whose lenght is equal to the levels of the factorial predictor, i. Each cell of the list is equal to the output of \code{delta.psy.probit} applied to a multivariable model whose baseline is level i of the factorial predictor.
}
\references{
Moscatelli A, Mezzetti M, Lacquaniti F (2012). Modelling Psychophysical Data at the Population-Level: The Generalized Linear Mixed Model. Journal of Vision, 12(11):26, 1-17.
}
\author{
Alessandro Moscatelli
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[lme4]{glmr}} for Generalized Linear Mixed Models (including random effects).\code{\link{delta.psy.probit}} for a univariable model. \code{\link{pseMer}} provides the bootstrap-based confidence intervals.
}
\examples{
#In two steps: Simulate a dataset with a categorical variable ("condition")
datafr.1 <- MERsimulate(fixeff = c(-7.5, 0.0875), nsubject = 6, constant = T)                               
levels(datafr.1$Subject) = c("S1", "S2", "S3", "S4", "S5", "S6")
datafr.1$condition = rep("A", 54)

datafr.2 <- MERsimulate(fixeff = c(-7, 0.0875),nsubject = 6, constant = T)                              
levels(datafr.2$Subject) = c("S1", "S2", "S3", "S4", "S5", "S6")
datafr.2$condition = rep("B", 54)

datafr = merge(datafr.1, datafr.2, all = T)
datafr$condition = as.factor(datafr$condition)

#How to estimate the PSE in the two condtions?
#1)use MERtreatment (note, you have to provide also the dataframe to the function)
formula.mod = cbind(Longer, Total - Longer) ~ X * condition + (1 + X| Subject)
mod1 <- glmer(formula = formula.mod, family = binomial(link = "probit"), data = datafr)
xplode.mod1 = xplode.mer(model = mod1, name.cont = "X", name.factor = "condition")
MERtreatment(xplode.mod1, datafr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Delta Method}
\keyword{Multivariable GLMM}% __ONLY ONE__ keyword per line
