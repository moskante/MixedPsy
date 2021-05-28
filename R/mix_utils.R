#' PSE/JND for univariable GLMM with Delta Method
#'
#' Estimate the Point of Subjective Equivalence (PSE), the Just Noticeable
#' Difference (JND) and the related Standard Errors for a univariate distribution 
#' (i.e. only one continuous predictor) by means of delta method.
#' 
#' @details \code{MixDelta} estimates PSE and JND of a univariable psychometric
#' function given an object of class \code{\link[lme4]{merMod}}.
#' The method only applies to GLMMs with only one continuous predictor and a 
#' \emph{probit} link function. Use \code{\link{MixTreatment}} for multivariable GLMMs with a factorial 
#' and a continuous predictor. \code{MixDelta} assumes that the first model coefficient is the intercept
#' and the second is the slope. The JND estimate assumes a \emph{probit} link function.
#'
#' @param xplode.obj an object of class \code{\link{xplode.obj}}.
#' @param alpha significance level of the confidence intervals. Default is 0.05. 
#'
#' @return \code{MixDelta} returns a list of length 1 including Estimate, Standard Error,
#' Inferior and Superior Confidence Interval of PSE and JND. Confidence Intervals
#' are computed as: \eqn{Estimate +/- z(1-(\alpha/2)) * Std.Error}.
#'
#' @note The delta method is based on the assumption of asymptotic normal distribution of the parameters estimates. 
#' This may result in an incorrect estimation. For a more reliable (but more time-consuming) 
#' bootstrap-based estimate, use \code{\link{pseMer}}.
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
#'  \code{\link{MixTreatment}} for a generalization of the parameter estimation for multivariable GLMM. 
#'  \code{\link{pseMer}} for bootstrap-based confidence intervals. 
#'  \code{\link{xplode}} for objects of class \code{\link{xplode.obj}}. 
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
    
    # copy all the variables in temporary objects
    pse <- -(xplode.obj$psychometrics[[1]]$intercept[1]/xplode.obj$psychometrics[[1]]$slope[1])
    slope <- xplode.obj$psychometrics[[1]]$slope[1]
    
    var.intercept <- xplode.obj$psychometrics[[1]]$intercept[2]
    var.slope <- xplode.obj$psychometrics[[1]]$slope[2]
    
    # cov(alpha, slope): for all pfs, is approximated to the cov(alpha1, slope1)
    cov.intercept.slope <- xplode.obj$psychometrics$pf1$cov
    
    # compute all the other variables
    var.pse <- (1/slope^2) * (var.intercept + (2 * pse * cov.intercept.slope) + (pse^2 * var.slope))  #PSE
    inferior.pse <- pse - (qnorm(1 - (alpha/2)) * sqrt(var.pse))
    superior.pse <- pse + (qnorm(1 - (alpha/2)) * sqrt(var.pse))
    
    jnd <- qnorm(0.75) * (1/slope)
    var.jnd <- (qnorm(0.75) * (-1/slope^2))^2 * var.slope  #JND
    inferior.jnd <- jnd - (qnorm(1 - (alpha/2)) * sqrt(var.jnd))
    superior.jnd <- jnd + (qnorm(1 - (alpha/2)) * sqrt(var.jnd))
    
    output <- matrix(rbind(c(pse, sqrt(var.pse), inferior.pse, superior.pse), 
                           c(jnd,sqrt(var.jnd), inferior.jnd, superior.jnd)), nrow = 2, 
                     dimnames = list(param <- c("pse","jnd"), statistics = c("Estimate", "Std. Error", "Inferior", "Superior")))
  }
  
  return(output)
}

#' PSE/JND for Multivariable GLMM with Delta Methods
#'
#' Estimate the Point of Subjective Equivalence (PSE), the Just Noticeable
#' Difference (JND) and the related Standard Errors for a multivariate distribution
#' (i.e. one continuous and one factorial predictor) by means of 
#' delta method.
#' The method applies to multivariable GLMM having a \emph{probit} link function.
#' The function is based on a recursive use of \code{\link[lme4]{glmer}} and
#' \code{\link{MixDelta}}
#'
#' @param xplode.obj an object of class \code{xplode.obj}. The fitted model
#' (object of class \code{\link[lme4]{merMod}}) from \code{xplode.obj} includes
#' one continuous predictor and one factorial predictor.
#' @param alpha significance level of the confidence intervals. Default is 0.05. 
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
#' @seealso \code{\link[lme4]{glmer}} for Generalized Linear Mixed Models. 
#' \code{\link{MixDelta}} for PSE and JND estimation from a univariable GLMM using delta method. 
#' \code{\link{pseMer}} for bootstrap-based confidence intervals of psychometric parameters. 
#'
#' @keywords DeltaMethod Multivariable GLMM
#'
#' @examples
#' library(lme4)
#' data(vibro_exp3)
#' formula.mod <- cbind(faster, slower) ~ speed * vibration + (1 + speed| subject)
#' mod <- glmer(formula = formula.mod, family = binomial(link = "probit"), data = vibro_exp3)
#' xplode.mod <- xplode(model = mod, name.cont = "speed", name.factor = "vibration")
#' MixTreatment(xplode.mod)
#'
#' @importFrom stats binomial contrasts<- contr.treatment
#' @export
#'
MixTreatment <- function(xplode.obj, alpha = 0.05) {
  
  datafr = xplode.obj$model.frame
  resp = split(datafr[[1]], col(datafr[[1]]))
  names(resp) = xplode.obj$response.colnames
  datafr = cbind(datafr, data.frame(resp))
  
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
    delta.par[[i]] = MixDelta(temp.xplode[[i]], alpha = alpha)
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

#' PSE/JND for GLMM Using Bootstrap Methods
#'
#' Estimates the Point of Subjective Equivalence (PSE), the Just Noticeable
#' Difference (JND) and the related Standard Errors by means of Bootstrap Method.
#' 
#' @param mer.obj An object of class \code{"\linkS4class{merMod}"}.
#' @param B integer: the number of bootstrap samples.
#' @param FUN An optional, custom made function to specify the required parameters to be estimated.
#' if NULL, \code{pseMer()} will estimate the PSE and the JND of a univariable GLMM.
#' @param alpha Significance level of the confidence interval.
#' @param ci.type A vector of character strings representing the type of intervals required. The value 
#' should be any subset of the values c("norm","basic", "stud", "perc", "bca") or simply "all" which will 
#' compute all five types of intervals. "perc" should be always included for the summary table.
#' @param beep Logical. If TRUE, a "ping" sound alerts that the simulation is complete.
#'
#' @return \code{pseMer} returns a list of length 3 including a summary table (Estimate, Standard Error,
#' Inferior and Superior Confidence Interval of the parameters) and the output of  \code{\link[lme4]{bootMer}} 
#' and \code{\link[boot]{boot.ci}} functions, for further analises. Confidence Intervals in the summary table are
#' based on the percentile method.
#'
#' @details \code{pseMer} estimates PSE and JND (and additional user defined paremters) from a 
#' fitted GLMM model (class \code{"\linkS4class{merMod}"}). 
#' The "ping" sound is provided by \code{\link[beepr]{beep}} function from the \code{beepr} package.
#' 
#' @note A first custom function was written in 2012 for the non-CRAN package MERpsychophisics,
#' based on the algorithm in Moscatelli et al. (2012). The current function is a simple wrapper
#' of \code{lme4::bootMer()} and \code{boot::boot.ci()} functions.
#' 
#' Increasing the nuber of bootstrap samples (\code{B}) makes the estimate more reliable. 
#' However, this will also increase the duration of the computation.
#'
#' @references
#' Moscatelli, A., Mezzetti, M., & Lacquaniti, F. (2012). Modeling psychophysical data 
#' at the population-level: The generalized linear mixed model. 
#' Journal of Vision, 12(11):26, 1-17. https://doi.org/10.1167/12.11.26
#' 
#' Bates, D., Mächler, M., Bolker, B., & Walker, S. (2015). Fitting Linear Mixed-Effects 
#' Models Using lme4. Journal of Statistical Software, 67(1), 51. https://doi.org/10.18637/jss.v067.i01
#'
#' @seealso
#' \code{\link[lme4]{bootMer}} from \code{lme4} package and \code{\link[boot]{boot.ci}} from \code{boot} package. 
#' 
#' @keywords Univariable Multivariable GLMM Bootstrap
#'
#' @examples
#' ## Example 1: estimate pse/jnd of a univariable GLMM
#' library(lme4)
#' data(vibro_exp3)
#' formula.mod1 <- cbind(faster, slower) ~ speed + (1 + speed| subject)
#' mod1 <- glmer(formula = formula.mod1, family = binomial(link = "probit"), 
#'               data = vibro_exp3[vibro_exp3$vibration == 0,])
#' \dontshow{BootEstim.1a <- pseMer(mod1, B = 5, ci.type = c("perc"))}
#' \donttest{BootEstim.1 <- pseMer(mod1, B = 100, ci.type = c("perc"))}
#' 
#' ## Example 2: specify custom parameters for bootstrap estimation of a 
#' # multivariate model
#' 
#' formula.mod2 <- cbind(faster, slower) ~ speed * vibration + (1 + speed| subject)
#' mod2 <- glmer(formula = formula.mod2, family = binomial(link = "probit"), 
#'                data = vibro_exp3)
#'               
#' fun2mod = function(mer.obj){
#' #allocate space: 4 parameters (jnd_0Hz, jnd_32Hz, pse_0Hz, pse_32Hz) j
#' jndpse = vector(mode = "numeric", length = 4)
#' names(jndpse) = c("jnd_0Hz","jnd_32Hz", "pse_0Hz", "pse_32Hz")
#' jndpse[1] = qnorm(0.75)/fixef(mer.obj)[2] #jnd_0Hz
#' jndpse[2] = qnorm(0.75)/(fixef(mer.obj)[2] + fixef(mer.obj)[4]) #jnd_32Hz
#' jndpse[3] = -fixef(mer.obj)[1]/fixef(mer.obj)[2] #pse_0Hz
#' jndpse[4] = -(fixef(mer.obj)[1] + fixef(mer.obj)[3])/(fixef(mer.obj)[2] 
#'                + fixef(mer.obj)[4]) #pse_32Hz
#' return(jndpse)
#' }
#' 
#' \donttest{BootEstim.2 = pseMer(mod2, B = 100, FUN = fun2mod)}
#' 
#' @export
#' @importFrom lme4 bootMer
#' @importFrom Matrix nearPD
#' @importFrom boot boot.ci
#' @importFrom beepr beep

pseMer <- function(mer.obj, B = 200, FUN = NULL, alpha = 0.05, 
                   ci.type = c("norm", "basic", "perc"), beep = F) {
  
  if (is.null(FUN)) {
    myfun <- function(mer.obj) {
      jndpse = vector(mode = "numeric", length = 2)
      names(jndpse) = c("JND", "PSE")
      jndpse[1] = qnorm(0.75)/fixef(mer.obj)[2]
      jndpse[2] = -fixef(mer.obj)[1]/fixef(mer.obj)[2]
      return(jndpse)
    }
  } else {
    myfun = match.fun(FUN)
  }
  
  np = length(myfun(mer.obj))
  parname = names(myfun(mer.obj))
  if (is.null(parname)) {
    print("Warning messages: Parameters have no names")
  }
  summary = matrix(NA, nrow = np, ncol = 3,
                   dimnames = list(parname, c("Estimate", "Inferior","Superior")))
  
  boot.samp <- lme4::bootMer(mer.obj, myfun, nsim = B)
  summary[, 1] = boot.samp$t0
  
  jndpseconf = vector(mode = "list", length = np)
  my.conf = 1 - alpha
  
  for (i in 1:np) {
    jndpseconf[[i]] <- boot::boot.ci(boot.samp, conf = my.conf, type = ci.type, index = i)
    if("perc" %in% ci.type){
      print(paste(parname[i], " 95% CI:", jndpseconf[[i]]$percent[4], "  ", jndpseconf[[i]]$percent[5]))
      summary[i, 2] = jndpseconf[[i]]$percent[4]
      summary[i, 3] = jndpseconf[[i]]$percent[5]
    }
  }
  
  if (beep == T) {
    beepr::beep()
  }
  out = list(summary, boot.samp, jndpseconf)
  return(out)
}

