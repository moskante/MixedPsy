#' PSE/JND from GLM Using Delta Method
#'
#' Estimate the Point of Subjective Equivalence (PSE), the Just Noticeable
#' Difference (JND) and the related Standard Errors by means of Delta Method.
#' The method only applies to GLMs (psychometric functions) with one continuous 
#' predictor and a \emph{probit} link function.
#'
#'
#' @param model the fitted psychometric function. An object of class \code{"glm"}.
#' @param alpha significance level of the confidence interval.Default is 0.05.
#'
#' @details \code{PsychDelta} estimates PSE and JND of a psychometric
#' function (object of class \code{"glm"}).
#' 
#' @return \code{PsychDelta} returns a matrix including Estimate, Standard Error,
#' Inferior and Superior Confidence Interval of PSE and JND. Confidence Intervals
#' are computed as: \eqn{Estimate +/- z(1-(\alpha/2)) * Std.Error}.
#'
#' @note The function assumes that the first model coefficient is the intercept
#' and the second is the slope. The estimate of the JND assumes a \emph{probit}
#' link function.
#'
#' @references
#' Faraggi, D., Izikson, P., & Reiser, B. (2003). Confidence intervals for the 50 per cent 
#' response dose. Statistics in medicine, 22(12), 1977-1988. https://doi.org/10.1002/sim.1368
#'
#' Moscatelli, A., Mezzetti, M., & Lacquaniti, F. (2012). Modeling psychophysical data 
#' at the population-level: The generalized linear mixed model. 
#' Journal of Vision, 12(11):26, 1-17. https://doi.org/10.1167/12.11.26
#'
#' @seealso \code{\link[stats]{glm}} for for Generalized Linear Models (without
#' random effects) and \code{\link[lme4]{glmer}} for Generalized Linear Mixed
#' Models (including random effects). \code{MixDelta} for GLMMs, (object of class
#' \code{\linkS4class{merMod}}). \code{\link{pseMer}} for bootstrap-based
#' confidence intervals.
#' 
#' @keywords Univariable GLM DeltaMethod
#'
#' @examples
#' #load 1 participant from simulated data
#' data.S1 <- subset(simul_data, Subject == "S1")
#' #fit a glm (probit link)
#' model.glm = glm(formula = cbind(Longer, Total - Longer) ~ X,
#' family = binomial(link = "probit"), data = data.S1)
#' PsychDelta(model.glm)
#' 
#' @importFrom stats vcov
#' @export
#'
PsychDelta <- function(model, alpha = 0.05) {

    pse <- -model$coef[1]/model$coef[2]
    BETA <- model$coef[2]

    var.alpha <- vcov(model)[1, 1]
    var.beta <- vcov(model)[2, 2]
    cov.alpha.beta <- vcov(model)[2, 1]

    var.pse <- (1/BETA^2) * (var.alpha + (2 * pse * cov.alpha.beta) + (pse^2 * var.beta))  #PSE
    inferior.pse <- pse - (qnorm(1 - (alpha/2)) * sqrt(var.pse))
    superior.pse <- pse + (qnorm(1 - (alpha/2)) * sqrt(var.pse))

    jnd <- qnorm(0.75) * (1/BETA)
    var.jnd <- (qnorm(0.75) * (-1/BETA^2))^2 * var.beta  #JND
    inferior.jnd <- jnd - (qnorm(1 - (alpha/2)) * sqrt(var.jnd))
    superior.jnd <- jnd + (qnorm(1 - (alpha/2)) * sqrt(var.jnd))

    output <- matrix(rbind(c(pse, sqrt(var.pse), inferior.pse, superior.pse), c(jnd, sqrt(var.jnd),
        inferior.jnd, superior.jnd)), nrow = 2, dimnames = list(param <- c("pse", "jnd"), statistics <- c("Estimate",
        "Std. Error", "Inferior", "Superior")))

    return(output)
}


#' Fitting Psychometric Functions
#'
#' Fit psychometric functions using either
#' \code{glm()} or \code{brglm()}; estimate PSE, JND, and the related
#' confidence intervals with Delta Method.
#'
#' @param ps.formula an object of class \code{\link[stats]{formula}}, such as \code{cbind(yes, no) ~ X}
#' @param ps.link link function for the binomial family of error distribution. Default is \code{"probit"}.
#' @param ps.data a data frame including the variables in the model.
#' @param br  logical. If TRUE, \code{\link[brglm]{brglm}} for bias reduction is used if values are equal to 0 or 1.
#'
#' @details Estimates are computed only for GLM of the type \code{F(Y) ~ X}, where X is a continuous
#' predictor. Std. Errors and 95\% confidence intervals
#' of PSE and JND are estimated via Delta Methods (Faraggi et al., 2003). Currently only working with \emph{probit} 
#' link function.
#'
#' @return \code{\link{PsychFunction}} returns a list including the fitted model,
#' the estimate of PSE and JND and a flag to indicate if \code{\link[brglm]{brglm}} was called.
#'
#' @references
#' Faraggi, D., Izikson, P., & Reiser, B. (2003). Confidence intervals for the 50 per cent 
#' response dose. Statistics in medicine, 22(12), 1977-1988. https://doi.org/10.1002/sim.1368
#'
#' Moscatelli, A., Mezzetti, M., & Lacquaniti, F. (2012). Modeling psychophysical data 
#' at the population-level: The generalized linear mixed model. 
#' Journal of Vision, 12(11):26, 1-17. https://doi.org/10.1167/12.11.26
#' 
#' @seealso \code{\link[stats]{glm}} for Generalized Linear Models; 
#' \code{\link[brglm]{brglm}} for fitting a GLM using bias reduction;
#' \code{\link{PsychPlot}} for plotting a psychometric function given a \code{\link[stats]{glm}} (or \code{\link[brglm]{brglm}}) object;
#' \code{\link{PsychShape}} for plotting a psychometric function given its PSE and JND.
#' 
#' @keywords Psychometric GLM DeltaMethod
#'
#' @examples
#' # simulate data from a single participant
#' data.S1 <- PsySimulate(fixeff = c(-7.5, 0.0875), nsubject = 1, constant = TRUE)
#' fit.S1 = PsychFunction(ps.formula = cbind(Longer, Total - Longer) ~ X,
#'                         ps.link = "probit", ps.data = data.S1)
#'                         
#' @importFrom brglm brglm
#' @importFrom stats glm predict terms
#' @export
#'
PsychFunction <-  function (ps.formula, ps.link, ps.data, br = F) {
  
  myfit = list()
  ps.terms = terms(ps.formula)
  model.glm = glm(formula = ps.formula, family = binomial(link = ps.link), 
                  data = ps.data)
  
  eps = 1e-15
  brflag = ifelse(1 - max(model.glm$fitted.values) <= eps & 
                    trunc(min(model.glm$fitted.values)) == 0, T, F)
  
  if (br == T & brflag == T) {
    model.glm = brglm(formula = ps.formula, family = binomial(link = ps.link), 
                      data = ps.data)
    warning("Binomial-response GLMs fitted using the bias-reduction method (brglm)")
  }
  
  myfit$model = model.glm
  
  if (ps.link == "probit") {
    myfit$estimate = PsychDelta(model.glm)
  } else {
    myfit$estimate = NA
    warning("Use the probit link function to get the estimate of the PSE and the JND")
  }
  
  myfit$info <- list(brflag = brflag)
  
  return(myfit)
}

#' Plotting Psychometric Functions given PSE and JND
#'
#' Plot a psychometric function with known PSE and JND
#' on an existing plot.
#'
#' @param pse,jnd point of subjective equivalende (PSE) and just noticeable difference (JND) of the desired psychometric function.
#' @param x.range a vector of length two specifying the range of the function.
#' @param ps.link a link function for the binomial family of error distribution (see Details).
#' @ps.type
#' @ps.size
#' @ps.color
#'
#' @details \code{PsychShape()} can be used to visualize the predicted results of a
#' psychophysical experiment or to plot a fitted psychometric function whose
#' values of pse and jnd are known. Currently only working with probit and logit
#' link function.
#'
#' @references
#' Moscatelli, A., Mezzetti, M., & Lacquaniti, F. (2012). Modeling psychophysical data 
#' at the population-level: The generalized linear mixed model. 
#' Journal of Vision, 12(11):26, 1-17. https://doi.org/10.1167/12.11.26
#' 
#' Knoblauch, K., & Maloney, L. T. (2012). Modeling psychophysical data in R (Vol. 32). 
#' Springer Science & Business Media.
#'
#' @seealso \code{\link[stats]{glm}} for for Generalized Linear Models.
#' \code{\link{PsychFunction}} for estimation of PSE and JND.
#'
#' @examples
#' p <- PsychShape(pse = 0, jnd = 6, x.range = c(-40, 40), ps.color = "gray", ps.size = 3)
#' p1 <- PsychShape(pse = 6, jnd = 6, x.range = c(-40, 40), ps.col = "black", addTo = p)
#' p2 <- PsychShape(pse = 6, jnd = 6, x.range = c(-40, 40), ps.col = "red", ps.link = "logit", ps.type = "dashed", addTo = NULL)
#'
#' @importFrom stats plogis
#' @import ggplot2
#' @export
#' 
PsychShape <- function(pse = 0, jnd = 1, x.range = c(NA, NA), ps.link = c("probit", "logit"), 
                       ps.type = "solid", ps.size = 1, ps.color = "black", addTo = NULL) {
  if(is.null(addTo)){
    p <- ggplot()
  }else{ 
    p <- addTo}
  
  x = pretty(x.range, 100)
  if (ps.link == "probit") {
    slope = qnorm(0.75) * (1/jnd)
    y <- pnorm(x, mean = pse, sd = 1/slope)
    
  } else if (ps.link == "logit") {
    slope = log(3) * (1/jnd)
    y = plogis(x, location = pse, scale = 1/slope)
  } else {
    warning("The function only works with probit and logit link function")
  }
  data <- data.frame(x,y)
  plot <- list()
  plot$ps <- geom_line(inherit.aes = FALSE, data = data, aes_string("x", "y"), 
                       linetype = ps.type, size = ps.size, color = ps.color)  
  
  print(p + plot)
}

#' Plotting Psychometric Functions given GLM
#'
#' Plot a psychometric function given an object of class \code{\link[stats]{glm}} or \code{\link[brglm]{brglm}}. 
#'
#' @param mod 
#' @param addTo 
#' @param showData 
#' @param ps.type line type.
#' @param ps.size line size.
#' @param ps.lab label.
#'
#' @details 
#'
#' @references
#' Moscatelli, A., Mezzetti, M., & Lacquaniti, F. (2012). Modeling psychophysical data 
#' at the population-level: The generalized linear mixed model. 
#' Journal of Vision, 12(11):26, 1-17. https://doi.org/10.1167/12.11.26
#' 
#' Knoblauch, K., & Maloney, L. T. (2012). Modeling psychophysical data in R (Vol. 32). 
#' Springer Science & Business Media.
#'
#' @seealso \code{\link[stats]{glm}} for for Generalized Linear Models.
#' \code{\link{PsychFunction}} for estimation of PSE and JND.
#'
#' @examples
#'
#' @import ggplot2
#' @export
#' 
PsychPlot <- function(mod, addTo = NULL, showData = TRUE,
                      ps.type = "solid", ps.size = 1, ps.lab = ""){
  if(is.null(addTo)){
    p <- ggplot()
  }else{ 
    p <- addTo}
  
  # GLM
  temp <- names(mod$model)
  xname <- temp[2]
  
  data <- mod$data
  data$y <- mod$y
  longData <- data.frame(pretty(data[[xname]], 1000))
  names(longData) <- xname
  longData$y <- predict(object = mod, newdata = longData, type = "response")
  
  plot <- list()
  
  plot$ps <- geom_line(inherit.aes = FALSE, data = longData, aes_string(xname, 'y', color = 'ps.lab'), 
                       linetype = ps.type, size = ps.size)  
  if(isTRUE(showData)){
    plot$data <- geom_point(inherit.aes = FALSE, data = data, aes_string(xname, 'y', color = 'ps.lab'))
  }
  
  if(mod$family[["link"]] == "probit"){
    psych <- PsychDelta(mod)
    plot$segment <- geom_segment(inherit.aes = FALSE, data = NULL, 
                                 aes(x = psych["pse", "Inferior"], xend = psych["pse", "Superior"], 
                                     y = 0.5, yend = 0.5, color = ps.lab),
                                 size = ps.size)
  } 
  print(p + plot)
  
}

