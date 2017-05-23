#' Extract values from a fitted GLMM object
#'
#' A function to extract values from an object of class \code{"\linkS4class{merMod}"}
#' (more specifically, from an object of subclass glmerMod).
#'
#'  @param model the GLMM fitted with \code{glmer}. An object of class
#'  \code{"\linkS4class{merMod}"}.
#'  @param name.cont  a string providing the name of the continuous predictor,
#'  as in the formula object of the fitted model
#'  @param name.factor a string providing the name of name of the categorical
#'  predictor, as in the formula object of the fitted model
#'  @param names.response optional. A string providing the name of name of the
#'  response variable, as in the formula object of the fitted model
#'  @param define.pf  optional. Specifies which parameter pertains to the intercept
#'  and which to the slope in the formula object.
#'
#'  @details For simplicity, several \code{MixedPsy} functions take as input an
#'  object of class \code{xplode} instead of an object of class
#'  \code{"\linkS4class{merMod}"}. Most of these functions assume by default that
#'  the continuous predictor is entered first in the formula object. It is possible
#'  to use a different order, this requires to specify which parameter pertains to
#'  the intercept and which to the slope, by changing \code{define.pf}.
#'
#'  @return An object of class \code{"\linkS4class{merMod}"} to be used with other
#'  \code{MixedPsy} functions.
#'
#'  @examples
#' datafr = PsySimulate(nsubjects = 10)
#' mod1 = glmer(formula = cbind(Longer, Total - Longer) ~ X + (1 | Subject),
#' family = binomial(link = "probit"), data = datafr)
#' xplode.mod1 = xplode(model = mod1, name.cont = "X")
#' MERdelta.probit(xplode.mod2)
#'
#' @export
xplode = function(model, name.cont = NA, name.factor = NA, names.response = NA, define.pf = list(pf1 = list(intercept = 1,
    slope = 2))) {

  #To DO: adjust at line 115 for binary
  # (single col) data
    xplode = vector("list", 27)
    names(xplode) = c("fixef", "fixef.vcov", "ranef", "n.psych.fun", "factor.col", "factor.colname","factor.levels", "factor.parnames",
        "cont.col", "cont.colname", "psychometrics", "ranef.stddev", "ranef.VarCov", "multirand",
        "Groups", "Gp", "Groups.levels", "Groups.colnames", "flist", "nobs", "size", "response.colnames",
        "model.frame", "model.matrix", "rps", "formula", "family")


    #fixed effects--------------------------
    xplode$fixef = fixef(model)

    # Variance-Covariance Matrix of the fixed effects
    xplode$fixef.vcov = vcov(model)
    
    #random effects
    xplode$ranef = ranef(model)
    
    # Number of psychometric functions
    xplode$n.psych.fun = length(define.pf)

    # Number of psychometric functions
    xplode$define.pf = define.pf

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
    xplode$psychometrics = vector(mode = "list", length = xplode$n.psych.fun)
    names.pf = character(length = xplode$n.psych.fun)

    for (i in 1:xplode$n.psych.fun) {
        names.pf[i] = paste("pf", i, sep = "")

        intercept.pointers = define.pf[[i]]$intercept
        slope.pointers = define.pf[[i]]$slope

        # estimate and Variance of the Intercept
        xplode$psychometrics[[i]]$intercept = numeric(length = 2)
        names(xplode$psychometrics[[i]]$intercept) = c("Estimate", "Variance")
        xplode$psychometrics[[i]]$intercept[1] = sum(fixef(model)[intercept.pointers])

        if (length(intercept.pointers) == 1) {
            xplode$psychometrics[[i]]$intercept[2] = xplode$fixef.vcov[intercept.pointers, intercept.pointers]
        } else {
            kombo.intercept = kombo(intercept.pointers)
            xplode$psychometrics[[i]]$intercept[2] = xplode$fixef.vcov[kombo.intercept[[1]][1, 1],
                kombo.intercept[[1]][2, 1]] + xplode$fixef.vcov[kombo.intercept[[1]][1, 2], kombo.intercept[[1]][2,
                2]] + (2 * xplode$fixef.vcov[kombo.intercept[[2]][1, 1], kombo.intercept[[2]][2,
                1]])
        }

        # estimate and Variance of the Slope
        xplode$psychometrics[[i]]$slope = numeric(length = 2)
        names(xplode$psychometrics[[i]]$slope) = c("Estimate", "Variance")
        xplode$psychometrics[[i]]$slope[1] = sum(fixef(model)[slope.pointers])
        if (length(slope.pointers) == 1) {
            xplode$psychometrics[[i]]$slope[2] = xplode$fixef.vcov[slope.pointers, slope.pointers]
        } else {
            kombo.slope = kombo(slope.pointers)
            xplode$psychometrics[[i]]$slope[2] = xplode$fixef.vcov[kombo.slope[[1]][1, 1], kombo.slope[[1]][2,
                1]] + xplode$fixef.vcov[kombo.slope[[1]][1, 2], kombo.slope[[1]][2, 2]] + (2 * xplode$fixef.vcov[kombo.slope[[2]][1,
                1], kombo.slope[[2]][2, 1]])
        }
    }
    names(xplode$psychometrics) = names.pf
    xplode$psychometrics$pf1$cov = xplode$fixef.vcov[name.cont, "(Intercept)"]

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

#FUNCTION: Draws a curve corresponding to a probit link function (beta > 0 and beta < 0)
CurveProbit = function(X, x.from, x.to){
  BETAplus = which(X[,2] > 0) 
  Xplus = X[BETAplus,]
  Xminus = X[-BETAplus,]
  
  if(nrow(Xplus) > 0){	
    apply(X = Xplus, MARGIN = 1,
          FUN = function(X) {curve(expr = pnorm(x, mean = -X[1]/X[2], sd = 1/X[2] ),
                                   from = x.from, to = x.to,
                                   col = X[3],
                                   add = T)}
    )}
  
  if(nrow(Xminus) > 0){
    apply(X = Xminus, MARGIN = 1,
          FUN = function(X) {curve(expr = 1 - pnorm(x, mean = -X[1]/X[2], sd = 1/abs(X[2] )),
                                   from = x.from, to = x.to,
                                   col = X[3],
                                   add = T)}
    )}
  return(BETAplus)
}