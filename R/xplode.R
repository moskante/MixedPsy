#' Extract values from a fitted GLMM object
#'
#' Extract values from an object of class \code{\linkS4class{merMod}}
#' (more specifically, from an object of subclass \code{glmerMod}). 
#'
#' @param model.obj The GLMM fitted with \code{glmer}. An object of class
#' \code{"\linkS4class{merMod}"}.
#' @param name.cont  A string providing the name of the continuous predictor,
#' as in the formula object of the fitted model
#' @param name.factor A string providing the name of name of the categorical
#' predictor, as in the formula object of the fitted model
#' @param names.response Optional. A string providing the name of name of the
#'  response variable, as in the formula object of the fitted model
#'  
#' @details For simplicity and maintenance reasons, several \code{MixedPsy} functions take as input an
#'  object of class \code{xplode} instead of an object of class \code{\linkS4class{merMod}}. 
#'
#'  
#' @keywords GLMM
#' 
#' @seealso \code{\link[lme4]{merMod-class}} and \code{\link[lme4]{glmer}}.
#' \code{\link{MixDelta}}, \code{\link{MixPlot}} for use of objects of class \code{xplode}.
#'
#' @examples
#' library(lme4)
#' multi.mod <- glmer(cbind(faster, slower) ~ speed * vibration  + (1 + speed| subject), 
#' family = binomial(link = "probit"), data = vibro_exp3)
#' xplode.mod <- xplode(multi.mod, name.cont = "speed", name.factor = "vibration")
#' MixPlot(xplode.mod)
#' MixDelta(xplode.mod)
#'
#' @import lme4
#' @export
#' 
xplode = function(model.obj, name.cont = NA, name.factor = NA, names.response = NA) {
  
  #To DO: adjust at line 115 for binary
  # (single col) data
  xplode = vector("list", 26)
  names(xplode) = c("fixef", "fixef.vcov", "ranef", "factor.col", "factor.colname","factor.levels", "factor.parnames",
                    "cont.col", "cont.colname", "psychometrics", "ranef.stddev", "ranef.VarCov", "multirand",
                    "Groups", "Gp", "Groups.levels", "Groups.colnames", "flist", "nobs", "size", "response.colnames",
                    "model.frame", "model.matrix", "rps", "formula", "family")
  
  
  #fixed effects--------------------------
  xplode$fixef = fixef(model.obj)
  
  # Variance-Covariance Matrix of the fixed effects
  xplode$fixef.vcov = vcov(model.obj)
  
  #random effects
  xplode$ranef = ranef(model.obj)
  
  # (the names of the parameters in the design matrix) to do: change from factor to treatment
  if (!is.na(name.factor)) {
    xplode$factor.col = which(names(model.frame(model.obj)) == name.factor)
    xplode$factor.colname = name.factor
    n.factors = length(levels(model.frame(model.obj)[, xplode$factor.col]))
    xplode$factor.levels = levels(model.frame(model.obj)[, xplode$factor.col])
    factor.parnames =  vector("character", n.factors)
    
    for (i in 1:n.factors) {
      xplode$factor.parnames[i] = paste(name.factor, levels(model.frame(model.obj)[, xplode$factor.col])[i],
                                        sep = "")
    }
  }
  # name of the continuous predictor
  xplode$cont.col = which(names(model.frame(model.obj)) == name.cont)
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
  xplode$psychometrics$intercept[1] = sum(fixef(model.obj)[intercept.pointers])
  
  xplode$psychometrics$intercept[2] = xplode$fixef.vcov[intercept.pointers, intercept.pointers]
  
  xplode$psychometrics$slope = numeric(length = 2)
  names(xplode$psychometrics$slope) = c("Estimate", "Variance")
  xplode$psychometrics$slope[1] = sum(fixef(model.obj)[slope.pointers])
  
  xplode$psychometrics$slope[2] = xplode$fixef.vcov[slope.pointers, slope.pointers]
  
  #names(xplode$psychometrics) = names.pf
  xplode$psychometrics$cov = xplode$fixef.vcov[name.cont, "(Intercept)"]
  
  # 2)random effects (univariate or multivariate)------------- Str. Dev.
  xplode$ranef.stddev <- as.numeric(attr(VarCorr(model.obj)[[1]], "stddev"))
  
  # Variance-Covariance Matrix of the random effects (two random effect and correlation)
  if (length(xplode$ranef.stddev) > 1) {
    xplode$ranef.VarCov = nearPD(VarCorr(model.obj)[[1]])$mat
    xplode$multirand = TRUE
  }else{
    xplode$multirand = FALSE
  }
  
  # 3)family xplode$family = family(model.obj) Bug fix 05.09.2014:
  xplode$family$family = summary(model.obj)$family
  xplode$family$link = summary(model.obj)$link
  
  # number of subjects and names
  xplode$Gp = getME(model.obj, name = "Gp")
  xplode$Groups.levels = levels(getME(model.obj, name = "flist")[[1]])
  xplode$flist = getME(model.obj, name = "flist")
  
  # find the column number in the frame specifiyng for the grouping factor
  xplode$Groups.colnames = names(xplode$flist)
  
  xplode$Groups = length(xplode$Groups.levels)
  
  xplode$nobs = nobs(model.obj)
  xplode$size = rowSums(model.response(model.frame(model.obj)))  #----------------------------> check this with binaries
  
  # 03.01.14: Check this
  if (is.na(names.response)) {
    xplode$response.colnames = colnames(model.response(model.frame(model.obj)))
  } else {
    xplode$response.colnames = names.response
  }
  
  xplode$model.frame = model.frame(model.obj)
  xplode$model.matrix = model.matrix(model.obj)
  
  # Rows per each Subject
  xplode$rps = table(getME(model.obj, "flist")[[1]])
  
  # formula
  xplode$formula = formula(model.obj)
  
  # modified on 14.11.2013 (yet not tested)
  class(xplode) = "xplode"
  
  # create the methds for the class, i.e., print, summary, ect. print.xplode <- function(x, ...){
  # cat('This is my vector:\n') cat(paste(x[1:5]), '...\n') }
  
  return(xplode)
}

