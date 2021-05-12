#' PSE/JND for Univariable GLMM Using Delta Methods
#'
#' Estimate the Point of Subjective Equivalence (PSE), the Just Noticeable
#' Difference (JND) and the related Standard Errors for an univariate distribution 
#' by means of Delta Method.
#' 
#' @details \code{MixDelta} estimates PSE and JND of a univariable psychometric
#' function (object of class \code{"glm"}).The method only applies to univariable GLMMs 
#'  having a \emph{probit} link function. Use \code{MixTreatment} for multivariable GLMMs.
#'
#' @param xplode.obj an object of class \code{xplode.obj} (univariable GLMMs).
#' @param alpha significance level of the confidence interval. Default is 0.05.
#'
#' @return \code{MixDelta} returns a list of length 1 including Estimate, Standard Error,
#' Inferior and Superior Confidence Interval of PSE and JND. Confidence Intervals
#' are computed as: \eqn{Estimate +/- z(1-(\alpha/2)) * Std.Error}.
#'
#' @note The function assumes that the first model coefficient is the intercept
#' and the second is the slope. The estimate of the JND assumes a \emph{probit}
#' link function.
#'
#' @references
#' Moscatelli, A., Mezzetti, M., & Lacquaniti, F. (2012). Modeling psychophysical data 
#' at the population-level: The generalized linear mixed model. 
#' Journal of Vision, 12(11):26, 1-17. https://doi.org/10.1167/12.11.26
#' 
#' Casella, G., & Berger, R. L. (2002). Statistical inference (2nd ed.). 
#' Pacific Grove, CA: Duxbury Press
#'
#' @seealso
#'  \code{\link{MixTreatment}} for univarible and multivariable GLMM. 
#'  \code{\link{pseMer}} for bootstrap-based confidence intervals. 
#'  \code{\link{xplode}} objects of class \code{xplode.obj}. 
#'
#' @keywords DeltaMethod Univariable GLMM
#' 
#' @examples
#' library(lme4)
#' data(vibro_exp3)
#' formula.mod <- cbind(faster, slower) ~ speed + (1 + speed| subject)
#' mod <- glmer(formula = formula.mod, family = binomial(link = "probit"), 
#'               data = vibro_exp3[vibro_exp3$vibration == 0,])
#' define.mod <- list(pf = list(intercept = 1, slope = 2))
#' xplode.mod <- xplode(model = mod, name.cont = "speed", define.pf = define.mod)
#' pse.jnd <- MixDelta(xplode.mod)
#' 
#' @importFrom stats qnorm
#' @importFrom grDevices palette
#' @export
#'
MixDelta <- function(xplode.obj, alpha = 0.05) {
  
  # check if link = probit
  if (xplode.obj$family$link != "probit") {
    output = NA
    print("Use a probit link function")
  } else {
    n.pf = length(xplode.obj$psychometrics)
    output = vector("list", length = n.pf)
    names(output) = names(xplode.obj$psychometrics)
    
    for (i in 1:n.pf) {
      # copy all the variables in temporary objects
      pse <- -(xplode.obj$psychometrics[[i]]$intercept[1]/xplode.obj$psychometrics[[i]]$slope[1])
      slope <- xplode.obj$psychometrics[[i]]$slope[1]
      
      var.intercept <- xplode.obj$psychometrics[[i]]$intercept[2]
      var.slope <- xplode.obj$psychometrics[[i]]$slope[2]
      
      # cov(alpha, slope): for all pfs, is approximated to the cov(alpha1, slope1)
      cov.intercept.slope <- xplode.obj$psychometrics$pf1$cov
      
      # compute all the other variables
      var.pse <- (1/slope^2) * (var.intercept + (2 * pse * cov.intercept.slope) + (pse^2 *
                                                                                     var.slope))  #PSE
      inferior.pse <- pse - (qnorm(1 - (alpha/2)) * sqrt(var.pse))
      superior.pse <- pse + (qnorm(1 - (alpha/2)) * sqrt(var.pse))
      
      jnd <- qnorm(0.75) * (1/slope)
      var.jnd <- (qnorm(0.75) * (-1/slope^2))^2 * var.slope  #JND
      inferior.jnd <- jnd - (qnorm(1 - (alpha/2)) * sqrt(var.jnd))
      superior.jnd <- jnd + (qnorm(1 - (alpha/2)) * sqrt(var.jnd))
      
      output[[i]] <- matrix(rbind(c(pse, sqrt(var.pse), inferior.pse, superior.pse), c(jnd,
                                                                                       sqrt(var.jnd), inferior.jnd, superior.jnd)), nrow = 2, dimnames = list(param <- c("pse",
                                                                                                                                                                         "jnd"), statistics = c("Estimate", "Std. Error", "Inferior", "Superior")))
    }
  }
  
  
  return(output)
}

#' PSE/JND for Multivariable GLMM Using Delta Methods
#'
#' Estimate the Point of Subjective Equivalence (PSE), the Just Noticeable
#' Difference (JND) and the related Standard Errors for a multivariate distribution by means of Delta Method.
#' The method applies to multivariable GLMM having a \emph{probit} link function.
#' The function is based on a recursive use of \code{glmer} and
#' \code{MixDelta}
#'
#' @param xplode.obj an object of class \code{xplode.obj}. The fitted model
#' (object of class \code{"\linkS4class{merMod}"}) from \code{xplode.obj} includes
#' one continuous predictor and one factorial predictor.
#' @param datafr  the data frame fitted with the GLMM model
#'
#' @details The function \code{MixTreatment} is based on a recursive use of
#' \code{glmer} and \code{PsychDelta} to multivariable GLMM including
#' continuous and factorial predictors. The same caveats of \code{PsychDelta}
#' apply (e.g., confidence interval based on normality assumption).
#'
#' @return A list, whose lenght is equal to the levels of the factorial predictor, i.
#' Each cell of the list is equal to the output of \code{delta.psy.probit} applied to
#' a multivariable model whose baseline is level i of the factorial predictor.
#'
#' @references
#' Moscatelli, A., Mezzetti, M., & Lacquaniti, F. (2012). Modeling psychophysical data 
#' at the population-level: The generalized linear mixed model. 
#' Journal of Vision, 12(11):26, 1-17. https://doi.org/10.1167/12.11.26
#'
#' @seealso \code{\link[lme4]{glmer}} for Generalized Linear Mixed Models (including
#' random effects).\code{\link{MixDelta}} for univariable model with delta method.
#' \code{\link{pseMer}} for bootstrap-based confidence intervals.
#'
#' @keywords DeltaMethod Multivariable GLMM
#'
#' @examples
#' library(lme4)
#' data(vibro_exp3)
#' formula.mod <- cbind(faster, slower) ~ speed * vibration + (1 + speed| subject)
#' mod <- glmer(formula = formula.mod, family = binomial(link = "probit"), data = vibro_exp3)
#' xplode.mod <- xplode(model = mod, name.cont = "speed", name.factor = "vibration")
#' MixTreatment(xplode.mod, vibro_exp3)
#'
#' @importFrom stats binomial contrasts<- contr.treatment
#' @export
#'
MixTreatment <- function(xplode.obj, datafr) {
  
  treat.lev = nlevels(xplode.obj$model.frame[, xplode.obj$factor.col])
  temp.models = delta.par = temp.xplode = vector("list", treat.lev)
  names(delta.par) = xplode.obj$factor.parnames
  
  for (i in 1:treat.lev) {
    contrasts(datafr[, which(names(datafr) == xplode.obj$factor.colname)]) = contr.treatment(treat.lev,
                                                                                             base = i)
    temp.models[[i]] = glmer(formula = xplode.obj$formula, family = binomial("probit"), data = datafr,
                             nAGQ = 1)
    temp.xplode[[i]] = xplode(temp.models[[i]], name.cont = xplode.obj$cont.colname, name.factor = xplode.obj$factor.colname,
                              define.pf = xplode.obj$define.pf)
    delta.par[[i]] = MixDelta(temp.xplode[[i]])[[1]]
  }
  return(delta.par)
}

#' Plotting univariable GLMM
#'
#' Plot binomial data and the fitted GLMM (object of class xplode).
#' 
#' @param xplode.obj an object of class xplode
#' @param pf integer: for multivariable GLMM including one factorial predictor,
#' the level number to be plotted 
#' @param p05line logical, should be an horizontal and a vertical line added? 
#' the horizontal line is fixed at P(Y = 1) = 0.5.
#' @param x.ref if p05line = T, this is the position of the vertical line on the x axis
#' @param x.range a vector of length two specifying the range for model predictions
#' @param col logical, if TRUE a different color will be used for different clusters/participants
#' @param x.label,y.label label for the x and the y axes. If not specified, x.labels = ""Stimulus Intensity", 
#' y.label = "Predicted Response"
#'   
#' @note The function is currently only working with GLMM including maximum three 
#' random effects (random intercept, random slope and covariance of the two)
#'  
#' @return a data.frame object including the intercept and slope for each participant
#' (algebraic sum of the fixed effects and the modes of the random effects) and the 
#' color number for the plot. 
#' 
#' @seealso \code{\link{xplode}} objects of class \code{xplode.obj}.
#' 
#' @keywords DeltaMethod GLMM Plotting
#'
#' @examples
#' library(lme4)
#' data(vibro_exp3)
#' formula.mod <- cbind(faster, slower) ~ speed + (1 + speed| subject)
#' mod <- glmer(formula = formula.mod, family = binomial(link = "probit"),
#'               data = vibro_exp3[vibro_exp3$vibration == 0,])
#' define.mod <- list(pf1 = list(intercept = 1, slope = 2))
#' xplode.mod <- xplode(model = mod, name.cont = "speed", define.pf = define.mod)
#' myplot <- MixPlot(xplode.mod, pf = 1,  p05line = FALSE, x.ref = 8.5, x.range = c(1,16),
#'                   col = TRUE, x.label = "Stimulus Speed", y.label = "Predicted Response")
#'
#' @export
#' 
MixPlot <- function(xplode.obj, pf = 1, p05line = F, x.range, x.ref,
                    col = F, x.label = "Stimulus Intensity", y.label = "Predicted Response"){
  
  #numebr of subjects
  nsubjects = xplode.obj$Groups
  
  #only one condition
  if(xplode.obj$n.psych.fun > 1){
    newdata = subset.data.frame(xplode.obj$model.frame,
                                xplode.obj$model.frame[ ,xplode.obj$factor.col] == xplode.obj$factor.levels[pf])
  }else{
    newdata = xplode.obj$model.frame
  }
  newdata$size = rowSums(newdata[,1])
  
  #colors and pch 
  if(col == T){
    if(length(palette()) < nsubjects){palette(rainbow(nsubjects))}		
    newdata$pointcolors = newdata[,xplode.obj$Groups.colnames]
    curvecolors = 1:nsubjects
  }else{
    newdata$pointcolors = 1
    curvecolors = rep(1, nsubjects)
  }
  
  #BLUPS (random modes)
  fixef.pf = cbind(xplode.obj$psychometrics[[pf]]$intercept["Estimate"],
                   xplode.obj$psychometrics[[pf]]$slope["Estimate"])
  
  if(xplode.obj$multirand == FALSE){
    estimates = data.frame(matrix(numeric(), nrow = nrow(xplode.obj$ranef[[1]]), ncol = 3))
    estimates[,1] = xplode.obj$ranef[[1]] + fixef.pf[1]
    estimates[,2] = rep(fixef.pf[2], times = nrow(estimates))
    estimates[,3] = curvecolors
  }else{
    estimates = t(apply(X = xplode.obj$ranef[[1]], MARGIN = 1, FUN = function(X){X + fixef.pf}))
    estimates = data.frame(cbind(estimates, curvecolors))
  }
  
  #see xplode$moddel.frame
  plot(y =  newdata[,1][,1]/newdata$size,
       x = newdata[,xplode.obj$cont.col],  xlim = x.range,
       xlab = x.label, ylim = c(0,1), ylab = y.label,
       col = newdata$pointcolors)		 													
  
  if (p05line == T){
    abline(h = 0.5, col = "lightgrey", lty = "dashed")
    abline(v = x.ref, col = "lightgrey", lty = "dashed")
  }
  
  CurveProbit(estimates, x.range[1], x.range[2])	
  
  palette("default")	
  return(estimates)
}