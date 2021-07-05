#' Extract values from a fitted GLMM object
#'
#' Extract values from an object of class \code{"\linkS4class{merMod}"}
#' (more specifically, from an object of subclass \code{glmerMod}).
#'
#' @param model The GLMM fitted with \code{glmer}. An object of class
#' \code{"\linkS4class{merMod}"}.
#' @param name.cont  A string providing the name of the continuous predictor,
#' as in the formula object of the fitted model
#' @param name.factor A string providing the name of name of the categorical
#' predictor, as in the formula object of the fitted model
#' @param names.response Optional. A string providing the name of name of the
#'  response variable, as in the formula object of the fitted model
#'  
#' @details For simplicity, several \code{MixedPsy} functions take as input an
#'  object of class \code{xplode} instead of an object of class
#'  \code{"\linkS4class{merMod}"}. Most of these functions assume by default that
#'  the continuous predictor is entered first in the formula object. It is possible
#'  to use a different order, this requires to specify which parameter pertains to
#'  the intercept and which to the slope, by changing \code{define.pf}.
#'
#' @return An object of class \code{\linkS4class{merMod}} to be used with other
#'  \code{MixedPsy} functions.
#'  
#' @keywords GLMM
#' 
#' @seealso \code{\link[lme4]{merMod-class}} and \code{\link[lme4]{glmer}} from package 
#' \code{\link[lme4]{lme4}} for objects of class ``\code{merMod}''.
#' \code{\link{MixDelta}}, \code{\link{MixDelta}} for use of objects of class \code{xplode.obj}.
#'
#' @examples
#' library(lme4)
#' datafr = PsySimulate(nsubjects = 10)
#' mod1 = glmer(formula = cbind(Longer, Total - Longer) ~ X + (1 | Subject),
#' family = binomial(link = "probit"), data = datafr)
#' xplode.mod1 = xplode(model = mod1, name.cont = "X")
#' MixDelta(xplode.mod1)
#'
#' @import lme4
#' @importFrom stats formula model.frame model.matrix model.response nobs
#' @export
#' 
xplode = function(model, name.cont = NA, name.factor = NA, names.response = NA) {
  
  #To DO: adjust at line 115 for binary
  # (single col) data
  xplode = vector("list", 26)
  names(xplode) = c("fixef", "fixef.vcov", "ranef", "factor.col", "factor.colname","factor.levels", "factor.parnames",
                    "cont.col", "cont.colname", "psychometrics", "ranef.stddev", "ranef.VarCov", "multirand",
                    "Groups", "Gp", "Groups.levels", "Groups.colnames", "flist", "nobs", "size", "response.colnames",
                    "model.frame", "model.matrix", "rps", "formula", "family")
  
  
  #fixed effects--------------------------
  xplode$fixef = fixef(model)
  
  # Variance-Covariance Matrix of the fixed effects
  xplode$fixef.vcov = vcov(model)
  
  #random effects
  xplode$ranef = ranef(model)
  
  # (the names of the parameters in the design matrix) to do: change from factor to treatment
  if (!is.na(name.factor)) {
    xplode$factor.col = which(names(model.frame(model)) == name.factor)
    xplode$factor.colname = name.factor
    n.factors = length(levels(model.frame(model)[, xplode$factor.col]))
    xplode$factor.levels = levels(model.frame(model)[, xplode$factor.col])
    factor.parnames =  vector("character", n.factors)
    
    for (i in 1:n.factors) {
      xplode$factor.parnames[i] = paste(name.factor, levels(model.frame(model)[, xplode$factor.col])[i],
                                        sep = "")
    }
  }
  # name of the continuous predictor
  xplode$cont.col = which(names(model.frame(model)) == name.cont)
  if (!is.na(xplode$cont.col)) {
    xplode$cont.colname = name.cont
  } else {
    print("Warning: The name of the continuous predictor is not among the names of the model.frame")
  }
  
  
  # re-arrange the fixed effects into n psychometric functions rename the n psychometric functions
  
  
  intercept.pointers = 1
  slope.pointers = 2
  
  # estimate and Variance of the Intercept
  xplode$psychometrics$intercept = numeric(length = 2)
  names(xplode$psychometrics$intercept) = c("Estimate", "Variance")
  xplode$psychometrics$intercept[1] = sum(fixef(model)[intercept.pointers])
  
  xplode$psychometrics$intercept[2] = xplode$fixef.vcov[intercept.pointers, intercept.pointers]
  
  xplode$psychometrics$slope = numeric(length = 2)
  names(xplode$psychometrics$slope) = c("Estimate", "Variance")
  xplode$psychometrics$slope[1] = sum(fixef(model)[slope.pointers])
  
  xplode$psychometrics$slope[2] = xplode$fixef.vcov[slope.pointers, slope.pointers]
  
  #names(xplode$psychometrics) = names.pf
  xplode$psychometrics$cov = xplode$fixef.vcov[name.cont, "(Intercept)"]
  
  # 2)random effects (univariate or multivariate)------------- Str. Dev.
  xplode$ranef.stddev <- as.numeric(attr(VarCorr(model)[[1]], "stddev"))
  
  # Variance-Covariance Matrix of the random effects (two random effect and correlation)
  if (length(xplode$ranef.stddev) > 1) {
    xplode$ranef.VarCov = nearPD(VarCorr(model)[[1]])$mat
    xplode$multirand = TRUE
  }else{
    xplode$multirand = FALSE
  }
  
  # 3)family xplode$family = family(model) Bug fix 05.09.2014:
  xplode$family$family = summary(model)$family
  xplode$family$link = summary(model)$link
  
  # number of subjects and names
  xplode$Gp = getME(model, name = "Gp")
  xplode$Groups.levels = levels(getME(model, name = "flist")[[1]])
  xplode$flist = getME(model, name = "flist")
  
  # find the column number in the frame specifiyng for the grouping factor
  xplode$Groups.colnames = names(xplode$flist)
  
  xplode$Groups = length(xplode$Groups.levels)
  
  xplode$nobs = nobs(model)
  xplode$size = rowSums(model.response(model.frame(model)))  #----------------------------> check this with binaries
  
  # 03.01.14: Check this
  if (is.na(names.response)) {
    xplode$response.colnames = colnames(model.response(model.frame(model)))
  } else {
    xplode$response.colnames = names.response
  }
  
  xplode$model.frame = model.frame(model)
  xplode$model.matrix = model.matrix(model)
  
  # Rows per each Subject
  xplode$rps = table(getME(model, "flist")[[1]])
  
  # formula
  xplode$formula = formula(model)
  
  # modified on 14.11.2013 (yet not tested)
  class(xplode) = "xplode"
  
  # create the methds for the class, i.e., print, summary, ect. print.xplode <- function(x, ...){
  # cat('This is my vector:\n') cat(paste(x[1:5]), '...\n') }
  
  return(xplode)
}

#' @importFrom utils combn
#FUNCTION: rearranges vector, Example: a = c(1,2,3,6) ka = kombo(a). used in xplode
kombo = function(vector) {
  n = length(vector)
  if (n > 1) {
    N = length(2:n)
    temp = vector("list", N)
    
    twice = matrix(NA, nrow = 2, ncol = n)
    for (i in 1:n) {
      twice[, i] = rep(vector[i], 2)
    }
    
    temp[[1]] = twice
    for (j in 2:n) {
      temp[[j]] = combn(vector, j)
    }
  } else {
    temp = vector
  }
  return(temp)
}