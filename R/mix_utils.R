# PSE/JND for univariable GLMM with Delta Method
#
# Estimate the Point of Subjective Equivalence (PSE), the Just Noticeable
# Difference (JND) and the related Standard Errors for a GLMM with
# only one continuous predictor (fixed intercept and slope), by means of delta method.
# 
# @details \code{MixDelta} estimates PSE and JND of a univariable psychometric
# function given an object of class \code{\link[lme4]{merMod}}.
# The method only applies to GLMMs with only one continuous predictor and a 
# \emph{probit} link function. Use \code{\link{MixDelta}} for multivariable GLMMs with a factorial 
# and a continuous predictor. \code{MixDelta} assumes that the first model coefficient is the intercept
# and the second is the slope. The JND estimate assumes a \emph{probit} link function.
#
# @param xplode.obj an object of class \code{\link{xplode.obj}}.
# @param alpha significance level of the confidence intervals. Default is 0.05. 
#
# @return \code{MixDelta} returns a list of length 1 including Estimate, Standard Error,
# Inferior and Superior Confidence Interval of PSE and JND. Confidence Intervals
# are computed as: \eqn{Estimate +/- z(1-(\alpha/2)) * Std.Error}.
#
# @note The delta method is based on the assumption of asymptotic normal distribution of the parameters estimates. 
# This may result in an incorrect estimation. For a more reliable (but more time-consuming) 
# bootstrap-based estimate, use \code{\link{pseMer}}.
#
# @references
# Moscatelli, A., Mezzetti, M., & Lacquaniti, F. (2012). Modeling psychophysical data 
# at the population-level: The generalized linear mixed model. 
# Journal of Vision, 12(11):26, 1-17. doi:10.1167/12.11.26
# 
# Casella, G., & Berger, R. L. (2002). Statistical inference (2nd ed.). 
# Pacific Grove, CA: Duxbury Press
#
# @seealso
#  \code{\link{MixDelta}} for a generalization of the parameter estimation for multivariable GLMM. 
#  \code{\link{pseMer}} for bootstrap-based confidence intervals. 
#  \code{\link{xplode}} for objects of class \code{\link{xplode.obj}}. 
#
# @keywords DeltaMethod Univariable GLMM
# 
# @examples
# library(lme4)
# data(vibro_exp3)
# formula.mod <- cbind(faster, slower) ~ speed + (1 + speed| subject)
# mod <- glmer(formula = formula.mod, family = binomial(link = "probit"), 
#               data = vibro_exp3[vibro_exp3$vibration == 0,])
# define.mod <- list(pf = list(intercept = 1, slope = 2))
# xplode.mod <- xplode(model = mod, name.cont = "speed", define.pf = define.mod)
# pse.jnd <- MixDelta(xplode.mod)
# 
#' @importFrom stats qnorm
#
MixFunction <- function(xplode.obj, alpha, p) {
  
  # check if link = probit
  if (xplode.obj$family$link != "probit") {
    output = NA
    print("Use a probit link function")
  } else {
    
    # copy all the variables in temporary objects
    pse <- -(xplode.obj$psychometrics$intercept[1]/xplode.obj$psychometrics$slope[1])
    slope <- xplode.obj$psychometrics$slope[1]
    
    var.intercept <- xplode.obj$psychometrics$intercept[2]
    var.slope <- xplode.obj$psychometrics$slope[2]
    
    # cov(alpha, slope): for all pfs, is approximated to the cov(alpha1, slope1)
    cov.intercept.slope <- xplode.obj$psychometrics$cov
    
    # compute all the other variables
    var.pse <- (1/slope^2) * (var.intercept + (2 * pse * cov.intercept.slope) + (pse^2 * var.slope))  #PSE
    inferior.pse <- pse - (qnorm(1 - (alpha/2)) * sqrt(var.pse))
    superior.pse <- pse + (qnorm(1 - (alpha/2)) * sqrt(var.pse))
    
    jnd <- qnorm(p) * (1/slope)
    var.jnd <- (qnorm(p) * (-1/slope^2))^2 * var.slope  #JND
    inferior.jnd <- jnd - (qnorm(1 - (alpha/2)) * sqrt(var.jnd))
    superior.jnd <- jnd + (qnorm(1 - (alpha/2)) * sqrt(var.jnd))
    
    output <- matrix(rbind(c(pse, sqrt(var.pse), inferior.pse, superior.pse), 
                           c(jnd,sqrt(var.jnd), inferior.jnd, superior.jnd)), nrow = 2, 
                     dimnames = list(param <- c("pse","jnd"), statistics = c("Estimate", "Std. Error", "Inferior", "Superior")))
  }
  
  return(output)
}

#' PSE/JND from GLMM Estimates using Delta Method
#'
#' Estimate Points of Subjective Equivalence (PSE), Just Noticeable
#' Differences (JND) and the related Standard Errors from a GLMM by means of 
#' delta method. 
#' The method applies to models with a \emph{probit} link function, one continuous predictor,
#' and one (optional) factorial predictor. 
#' 
#' @param xplode.obj an object of class \code{xplode.obj}. The fitted model
#' (object of class \code{\link[lme4]{merMod}}, specifically of subclass \code{glmerMod}) includes
#' one continuous predictor and one (optional) factorial predictor.
#' @param alpha significance level of the confidence intervals. Default is 0.05 (value for 95\% confidence interval).
#' @param p probability value relative to the JND upper limit. Default is 0.75 (value for 50\% JND).
#'
#' @details When the model includes a factorial predictor, the function is based on a recursive use of
#' \code{\link[lme4]{glmer}} and re-order of levels of the factorial predictor. 
#' The JND estimate assumes a \emph{probit} link function.
#'
#' @return A matrix including estimate, standard error, inferior and superior bounds of the confidence interval of PSE and JND. 
#' If a factorial predictor is included in the model, the function returns a list, each item containing a matrix 
#' for the estimates relative to a level of the predictor.
#' 
#' @note The delta method is based on the assumption of asymptotic normal distribution of the parameters estimates. 
#' This may result in an incorrect variance estimation. For a more reliable (but more time-consuming) estimation 
#' based on bootstrap method, use \code{\link{pseMer}}.
#'
#' @references
#' Moscatelli, A., Mezzetti, M., & Lacquaniti, F. (2012). Modeling psychophysical data 
#' at the population-level: The generalized linear mixed model. 
#' Journal of Vision, 12(11):26, 1-17. doi:10.1167/12.11.26
#'
#' Casella, G., & Berger, R. L. (2002). Statistical inference (2nd ed.). 
#' Pacific Grove, CA: Duxbury Press

#' @seealso \code{\link[lme4]{glmer}} for fitting Generalized Linear Mixed Models. 
#' \code{\link{xplode}} for interfacing values from a fitted GLMM to \code{MixedPsy} functions. 
#' \code{\link{pseMer}} for bootstrap-based confidence intervals of psychometric parameters. 
#'
#' @keywords DeltaMethod GLMM
#'
#' @examples
#' library(lme4)
#' 
#' #univariable GLMM (one continuous predictor)
#' mod.uni = glmer(formula = cbind(Longer, Total - Longer) ~ X + (1 | Subject), 
#' family = binomial(link = "probit"), data = simul_data)
#' xplode.uni = xplode(model = mod.uni, name.cont = "X")
#' MixDelta(xplode.uni)
#'
#' #multivariable GLMM (one continuous and one factorial predictor)
#' mod.multi <- glmer(cbind(faster, slower) ~ speed * vibration + (1 + speed| subject), 
#' family = binomial(link = "probit"), data = vibro_exp3)
#' xplode.multi <- xplode(model = mod.multi, name.cont = "speed", name.factor = "vibration")
#' MixDelta(xplode.multi)
#'
#' @importFrom stats update contrasts contr.treatment binomial
#' @export
#'
MixDelta <- function(xplode.obj, alpha = 0.05, p = 0.75) {
  
  datafr = xplode.obj$model.frame
  temp.formula = xplode.obj$formula
  
  if(is.matrix(datafr[[1]])){
    resp = split(datafr[[1]], col(datafr[[1]]))
    names(resp) = c("y1", "y2")
    datafr = cbind(datafr, data.frame(resp))
    temp.formula = update(temp.formula, cbind(y1,y2) ~ .)
  }
  
  treat.lev = length(xplode.obj$factor.levels)
  
  if (treat.lev > 0){
    temp.models = delta.par = temp.xplode = vector("list", treat.lev)
    names(delta.par) = xplode.obj$factor.parnames
    
    for (i in 1:treat.lev) {
      contrasts(datafr[, which(names(datafr) == xplode.obj$factor.colname)]) = contr.treatment(treat.lev,
                                                                                               base = i)
      temp.models[[i]] = glmer(formula = temp.formula, family = binomial("probit"), data = datafr,
                               nAGQ = 1)
      temp.xplode[[i]] = xplode(temp.models[[i]], name.cont = xplode.obj$cont.colname, name.factor = xplode.obj$factor.colname)
      delta.par[[i]] = MixFunction(temp.xplode[[i]], alpha = alpha, p = p)
    }
  }else{
    delta.par <- MixFunction(xplode.obj, alpha = alpha, p = p)
  }
  return(delta.par)
}

#' Plot Individual Responses from GLMM
#'
#' Plot response curve for each individual in a population sample, given a GLMM with one continuous predictor 
#' and one (optional) factorial predictor. If the factorial predictor is specified, the response is plotted separately for each
#' individual and each predictor level. 
#' 
#' @param xplode.obj an object of class \code{\link{xplode}}.
#' @param facet_by optional. A string specifying the name of the faceting variable (either the participant identification 
#' or the factorial predictor). 
#' @param showData logical, defines if proportion of binomial responses for each stimulus level are presented. Default is TRUE.
#' 
#' @details If the model includes only a continuous predictor, the figure consist of a single panel, and each individual's 
#' response is assigned a different color. If a factorial predictor is included in the model, the faceting variable can be either the 
#' participant identification or the factorial predictor. By default, each panel shows an individual's response, different levels of 
#' the factorial predictor are coded by color.
#' 
#' @return \code{MixPlot} returns a \code{\link[ggplot2]{ggplot}} object. 
#' 
#' @seealso \code{\link{xplode}} for objects of class \code{xplode}.
#' \code{\link[ggplot2]{ggplot2}} for creating data visualizations. 
#' \code{\link{PsychPlot}} for plotting a a psychometric function from a GLM. 
#' 
#' @examples
#' library(lme4)
#' mod.multi <- glmer(cbind(faster, slower) ~ speed * vibration + (1 + speed| subject), 
#' family = binomial(link = "probit"), data = vibro_exp3)
#' xplode.multi <- xplode(model = mod.multi, name.cont = "speed", name.factor = "vibration")
#' 
#' MixPlot(xplode.multi)
#' #alternative visualization
#' MixPlot(xplode.multi, facet_by = "vibration", showData = FALSE)
#'
#' @keywords GLMM plot
#' 
#' @import ggplot2
#' @importFrom stats update binomial predict
#' 
#' @export
#' 

MixPlot <- function(xplode.obj, facet_by= NULL, showData = TRUE){

  xname = xplode.obj$cont.colname
  yname = xplode.obj$response.colnames[1]
  
  factorname <- xplode.obj$factor.colname
  if(is.null(factorname)){
    factorlevels <- NA
  }else{
    factorlevels <- xplode.obj$factor.levels
  }
  
  groupname = xplode.obj$Groups.colnames
  grouplevels = xplode.obj$Groups.levels
  
  temp.formula = xplode.obj$formula
  
  temp.data = as.data.frame(xplode.obj$model.frame)
  if(is.matrix(temp.data[[1]])){
    resp = split(temp.data[[1]], col(temp.data[[1]]))
    names(resp) = c("y1", "y2")
    temp.data = cbind(temp.data, data.frame(resp))
    temp.formula = update(temp.formula, cbind(y1,y2) ~ .)
  }
  temp.data$size = xplode.obj$size
  
  temp.model = glmer(formula = temp.formula, family = binomial("probit"), data = temp.data,
                           nAGQ = 1)
  
  longData = expand.grid(cont = pretty(temp.data[[xname]], 1000),
                          group = grouplevels, 
                          factor = factorlevels)
  names(longData) = c(xname,groupname,factorname)
  longData$y = predict(object = temp.model, newdata = longData, type = "response")
  
  
  p = ggplot(longData, aes_string(xname,"y")) +
    theme_bw() + 
    theme(panel.grid = element_blank(), strip.background = element_blank())
  
  plot = list()
  if(is.null(factorname)){
    colorname = groupname
  }else if (is.null(facet_by) || facet_by ==  groupname){
    colorname = factorname
    plot$facet = facet_wrap(c(groupname))
  }else if(facet_by ==  factorname){
    colorname = groupname
    plot$facet = facet_wrap(c(factorname))
  }
  
  plot$line = geom_line(aes_string(color = colorname))
  if(isTRUE(showData)){
    plot$data = geom_point(data = temp.data, aes_string(xname, "y1/size", color = colorname))
  }
  
  print(p + plot)
  
}


#' PSE/JND from GLMM Estimates Using Bootstrap Method
#'
#' Estimates the Point of Subjective Equivalence (PSE), the Just Noticeable
#' Difference (JND) and the related Standard Errors by means of Bootstrap Method, 
#' given an object of class \code{\link[lme4]{merMod}}.
#' 
#' @param mer.obj an object of class \code{\linkS4class{merMod}}.
#' @param B integer. Number of bootstrap samples.
#' @param FUN an optional, custom made function to specify the required parameters to be estimated.
#' If NULL, \code{pseMer} estimates PSE and 50\%JND of a univariable GLMM with a single intercept and slope.
#' @param alpha significance level of the confidence intervals. Default is 0.05 (95\% confidence interval).
#' @param ci.type vector of character strings representing the type of intervals required. The value 
#' should be any subset of the values accepted by \code{\link[boot]{boot.ci}}: c("norm","basic", "stud", "perc", "bca"). 
#' Specify "all" for all five types of intervals. "perc" should be always included for the summary table.
#' @param beep logical. If TRUE, a "ping" sound alerts that the simulation is complete. Default is FALSE.
#'
#' @return \code{pseMer} returns a list of length 3 including a summary table (estimate,
#' inferior and superior bounds of the confidence interval), the output of \code{\link[lme4]{bootMer}}, and that of 
#' \code{\link[boot]{boot.ci}}, for further analyses. Confidence intervals in the summary table are
#' based on the percentile method.
#'
#' @details \code{pseMer} estimates PSE and JND (and additional user defined parameters) from a 
#' fitted GLMM model (class \code{\linkS4class{merMod}}). 
#' 
#' @note A first custom function was written in 2012 for the non-CRAN package MERpsychophisics,
#' based on the algorithm in Moscatelli et al. (2012). The current function is a wrapper
#' of function \code{\link[lme4]{bootMer}} and \code{\link[boot]{boot.ci}}. 
#' 
#' Increasing the number of bootstrap samples (\code{B}) makes the estimate more reliable. 
#' However, this will also increase the duration of the computation.
#'
#' @references
#' Moscatelli, A., Mezzetti, M., & Lacquaniti, F. (2012). Modeling psychophysical data 
#' at the population-level: The generalized linear mixed model. 
#' Journal of Vision, 12(11):26, 1-17. doi:10.1167/12.11.26
#' 
#' Bates, D., Mächler, M., Bolker, B., & Walker, S. (2015). Fitting Linear Mixed-Effects 
#' Models Using lme4. Journal of Statistical Software, 67(1), 51. https://doi.org/10.18637/jss.v067.i01
#'
#' @seealso
#' \code{\link[lme4]{bootMer}} and \code{\link[boot]{boot.ci}} for estimation of confidence intervals with the bootstrap method. 
#' \code{\link{MixDelta}} for confidence intervals with delta method. 
#' 
#' @keywords Bootstrap GLMM
#'
#' @examples
#' library(lme4)
#' #example 1: univariable GLMM
#' mod.uni = glmer(formula = cbind(Longer, Total - Longer) ~ X + (1 | Subject),
#' family = binomial(link = "probit"), data = simul_data)
#' \dontshow{BootEstim.1 <- pseMer(mod.uni, B = 5, ci.type = c("perc"))}
#' \donttest{BootEstim.uni <- pseMer(mod.uni, B = 100, ci.type = c("perc"))}
#' 
#' #example 2: specify custom parameters for multivariable model
#' mod.multi <- glmer(cbind(faster, slower) ~ speed * vibration + (1 + speed| subject), 
#' family = binomial(link = "probit"), data = vibro_exp3)
#'               
#' fun2mod = function(mer.obj){
#' #allocate space: 4 parameters (jnd_A, jnd_B, pse_A, pse_B)
#' jndpse = vector(mode = "numeric", length = 4)
#' names(jndpse) = c("pse_0", "pse_32","jnd_0", "jnd_32")
#' jndpse[1] = -fixef(mer.obj)[1]/fixef(mer.obj)[2] #pse_0
#' jndpse[2] = -(fixef(mer.obj)[1]+fixef(mer.obj)[3])/(fixef(mer.obj)[2]+ fixef(mer.obj)[4]) #pse_0
#' jndpse[3] = qnorm(0.75)/fixef(mer.obj)[2] #jnd_0
#' jndpse[4] = qnorm(0.75)/(fixef(mer.obj)[2]+ fixef(mer.obj)[4]) #jnd_32
#' return(jndpse)
#' }
#'  
#' \donttest{BootEstim.multi = pseMer(mod.multi, B = 100, FUN = fun2mod)}
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
  
  
  boot.samp <- bootMer(mer.obj, myfun, nsim = B, .progress = "txt")
  summary[, 1] = boot.samp$t0
  
  jndpseconf = vector(mode = "list", length = np)
  my.conf = 1 - alpha
  cat("\n")
  for (i in 1:np) {
    jndpseconf[[i]] <- boot.ci(boot.samp, conf = my.conf, type = ci.type, index = i)
    if("perc" %in% ci.type){
      print(paste(parname[i], " ", 100*(1-alpha),"% CI:", jndpseconf[[i]]$percent[4], "  ", jndpseconf[[i]]$percent[5]))
      summary[i, 2] = jndpseconf[[i]]$percent[4]
      summary[i, 3] = jndpseconf[[i]]$percent[5]
    }else{
      if(i == 1){
        print("For summary table, include percentile method")
      }
      
    }
  }
  
  if (beep == T) {
    beepr::beep()
  }
  out = list(summary = summary, boot.samp = boot.samp, CI = jndpseconf)
  return(out)
}

