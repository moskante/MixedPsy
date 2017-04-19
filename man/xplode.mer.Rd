\name{xplode.mer}
\alias{xplode.mer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xplode.mer(model, name.cont = NA, name.factor = NA, names.response = NA, define.pf = list(pf1 = list(intercept = 1, slope = 2)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{name.cont}{
%%     ~~Describe \code{name.cont} here~~
}
  \item{name.factor}{
%%     ~~Describe \code{name.factor} here~~
}
  \item{names.response}{
%%     ~~Describe \code{names.response} here~~
}
  \item{define.pf}{
%%     ~~Describe \code{define.pf} here~~
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
function (model, name.cont = NA, name.factor = NA, names.response = NA, 
    define.pf = list(pf1 = list(intercept = 1, slope = 2))) 
{
    xplode = vector("list", 25)
    names(xplode) = c("fixef", "fixef.vcov", "n.psych.fun", "factor.col", 
        "factor.colname", "factor.parnames", "cont.col", "cont.colname", 
        "psychometrics", "ranef.stddev", "ranef.VarCov", "multirand", 
        "Groups", "Gp", "Groups.levels", "Groups.colnames", "flist", 
        "nobs", "size", "response.colnames", "model.frame", "model.matrix", 
        "rps", "formula", "family")
    xplode$fixef = fixef(model)
    xplode$fixef.vcov = vcov(model)
    xplode$n.psych.fun = length(define.pf)
    if (!is.na(name.factor)) {
        xplode$factor.col = which(names(model.frame(model)) == 
            name.factor)
        xplode$factor.colname = name.factor
        n.factors = length(levels(model.frame(model)[, xplode$factor.col]))
        factor.parnames = vector("character", n.factors)
        for (i in 1:n.factors) {
            xplode$factor.parnames[i] = paste(name.factor, levels(model.frame(model)[, 
                xplode$factor.col])[i], sep = "")
        }
    }
    xplode$cont.col = which(names(model.frame(model)) == name.cont)
    if (!is.na(xplode$cont.col)) {
        xplode$cont.colname = name.cont
    }
    else {
        print("Warning: The name of the continuous predictor is not among the names of the model.frame")
    }
    xplode$psychometrics = vector(mode = "list", length = xplode$n.psych.fun)
    names.pf = character(length = xplode$n.psych.fun)
    for (i in 1:xplode$n.psych.fun) {
        names.pf[i] = paste("pf", i, sep = "")
        intercept.pointers = define.pf[[i]]$intercept
        slope.pointers = define.pf[[i]]$slope
        xplode$psychometrics[[i]]$intercept = numeric(length = 2)
        names(xplode$psychometrics[[i]]$intercept) = c("Estimate", 
            "Variance")
        xplode$psychometrics[[i]]$intercept[1] = sum(fixef(model)[intercept.pointers])
        if (length(intercept.pointers) == 1) {
            xplode$psychometrics[[i]]$intercept[2] = xplode$fixef.vcov[intercept.pointers, 
                intercept.pointers]
        }
        else {
            kombo.intercept = kombo(intercept.pointers)
            xplode$psychometrics[[i]]$intercept[2] = xplode$fixef.vcov[kombo.intercept[[1]][1, 
                1], kombo.intercept[[1]][2, 1]] + xplode$fixef.vcov[kombo.intercept[[1]][1, 
                2], kombo.intercept[[1]][2, 2]] + (2 * xplode$fixef.vcov[kombo.intercept[[2]][1, 
                1], kombo.intercept[[2]][2, 1]])
        }
        xplode$psychometrics[[i]]$slope = numeric(length = 2)
        names(xplode$psychometrics[[i]]$slope) = c("Estimate", 
            "Variance")
        xplode$psychometrics[[i]]$slope[1] = sum(fixef(model)[slope.pointers])
        if (length(slope.pointers) == 1) {
            xplode$psychometrics[[i]]$slope[2] = xplode$fixef.vcov[slope.pointers, 
                slope.pointers]
        }
        else {
            kombo.slope = kombo(slope.pointers)
            xplode$psychometrics[[i]]$slope[2] = xplode$fixef.vcov[kombo.slope[[1]][1, 
                1], kombo.slope[[1]][2, 1]] + xplode$fixef.vcov[kombo.slope[[1]][1, 
                2], kombo.slope[[1]][2, 2]] + (2 * xplode$fixef.vcov[kombo.slope[[2]][1, 
                1], kombo.slope[[2]][2, 1]])
        }
    }
    names(xplode$psychometrics) = names.pf
    xplode$psychometrics$pf1$cov = xplode$fixef.vcov[name.cont, 
        "(Intercept)"]
    xplode$ranef.stddev <- as.numeric(attr(VarCorr(model)[[1]], 
        "stddev"))
    if (length(xplode$ranef.stddev) > 1) {
        xplode$ranef.VarCov = nearPD(VarCorr(model)[[1]])$mat
        xplode$multirand = TRUE
    }
    xplode$family$family = summary(model)$family
    xplode$family$link = summary(model)$link
    xplode$Gp = getME(model, name = "Gp")
    xplode$Groups.levels = levels(getME(model, name = "flist")[[1]])
    xplode$flist = getME(model, name = "flist")
    xplode$Groups.colnames = names(xplode$flist)
    xplode$Groups = length(xplode$Groups.levels)
    xplode$nobs = nobs(model)
    xplode$size = rowSums(model.response(model.frame(model)))
    if (is.na(names.response)) {
        xplode$response.colnames = colnames(model.response(model.frame(model)))
    }
    else {
        xplode$response.colnames = names.response
    }
    xplode$model.frame = model.frame(model)
    xplode$model.matrix = model.matrix(model)
    xplode$rps = table(getME(model, "flist")[[1]])
    xplode$formula = formula(model)
    class(xplode) = "xplode"
    return(xplode)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
