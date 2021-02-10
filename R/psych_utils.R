#' PSE/JND for univariable GLM Using Delta Method
#'
#' Estimate the Point of Subjective Equivalence (PSE), the Just Noticeable
#' Difference (JND) and the related Standard Errors by means of Delta Method.
#' The method only applies to univariable GLMs (psychometric functions) having
#' a \emph{probit} link function.
#'
#'
#' @param model the fitted psychometric function. An object of class \code{"glm"}.
#' @param alpha significance level of the confidence interval.
#'
#' @details \code{PsychDelta} estimates PSE and JND of a univariable psychometric
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
#' Models (including random effects). \code{MixDelta} and \code{MixTreatment}
#' for univarible and multivariable GLMM, respectively (object of class
#' \code{"\linkS4class{merMod}"}). \code{\link{pseMer}} for bootstrap-based
#' confidence intervals.
#' 
#' @keywords Univariable GLM DeltaMethod
#'
#' @examples
#' #load simulated data
#' data(simul_data)
#' #fit a glm (probit link)
#' model.glm = glm(formula = cbind(Longer, Total - Longer) ~ X,
#' family = binomial(link = "probit"), data = simul_data)
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


#' Fitting and Plotting Psychometric Functions
#'
#' Fit psychometric functions using either
#' \code{glm()} or \code{brglm()}, estimate PSE, JND and the related
#' confidence intervals, and draw the curve on an existing plot.
#'
#' @param ps.formula an object of class ``formula'', such as \code{cbind(yes, no) ~ X}
#' @param ps.link a link function for the binomial family of error distribution.
#' See `Details
#' @param ps.data a data frame including the variables in the model
#' @param ps.lines logical. If TRUE, model predictions and confidence intervals of
#' the PSE will be added to an existing plot
#' @param ps.x optionally, a data frame in which to look for variables with which
#' to predict. See `Details'
#' @param x.range a vector of length two specifying the range for model predictions
#' @param ps.col color of the lines to be plotted
#' @param ps.lty line type
#' @param ps.lwd line width
#' @param br  logical. If TRUE, brglm is used if fitted values are equal to 0 or 1
#'
#' @details If \code{lines = TRUE}, the function draws model predictions on an existing plot.
#' Only for univariable glm of the type \code{F(Y) ~ X}, where X is a continuous
#' predictor. If \code{ps.x} is empty, the new data frame is a vector of length = 1000,
#' whose range is specified from \code{x.range}. Std. Errors and 95\% confidence intervals
#' of the PSE and JND are estimated via Delta Methods, see Faraggi et al. (2003).
#'
#' @return a list including the fitted glm (or \code{brglm}),
#' the estimate of PSE and JND and a flag to indicate if \code{brglm} was called.
#'
#' @references
#' Faraggi, D., Izikson, P., & Reiser, B. (2003). Confidence intervals for the 50 per cent 
#' response dose. Statistics in medicine, 22(12), 1977-1988. https://doi.org/10.1002/sim.1368
#'
#' Moscatelli, A., Mezzetti, M., & Lacquaniti, F. (2012). Modeling psychophysical data 
#' at the population-level: The generalized linear mixed model. 
#' Journal of Vision, 12(11):26, 1-17. https://doi.org/10.1167/12.11.26
#' 
#' @seealso \code{\link[stats]{glm}} for for Generalized Linear Models.
#' \code{\link{PsychShape}} for plotting psychometric function of given PSE and JND
#' 
#' @keywords Plotting Psychometric GLM DeltaMethod
#'
#' @examples
#' # simulate data from a single participant
#' datafr.S1 <- PsySimulate(fixeff = c(-7.5, 0.0875), nsubject = 1, constant = TRUE)
#' #fit a glm (probit link)
#' model.glm = glm(formula = cbind(Longer, Total - Longer) ~ X,
#' family = binomial(link = "probit"), data = datafr.S1)
#'
#' #fit psychometric function single-subject data and draw on existing plot
#' plot(Longer/Total ~ X, data = datafr.S1)
#' fit.S1 = PsychFunction(ps.formula = cbind(Longer, Total - Longer) ~ X,
#'                         ps.link = "probit", ps.data = datafr.S1,
#'                         x.range = c(40, 120), ps.lines = TRUE)
#'                         
#' @importFrom brglm brglm
#' @importFrom stats glm predict terms
#' @importFrom graphics lines segments
#' @export
#'
PsychFunction <- function(ps.formula, ps.link, ps.data, x.range = c(NA, NA), ps.x = NA, ps.lines = F,
                          ps.col = "black", ps.lty = "dashed", ps.lwd = 1, br = F) {
  myfit = vector("list", 3)
  ps.terms = terms(ps.formula)

  model.glm = glm(formula = ps.formula, family = binomial(link = ps.link), data = ps.data)

  eps = 1e-15
  brflag = ifelse(1-max(model.glm$fitted.values) <= eps & trunc(min(model.glm$fitted.values)) ==
                    0, T, F)
  myfit[[3]] = brflag

  if (br == T & brflag == T) {
    model.glm = brglm(formula = ps.formula, family = binomial(link = ps.link), data = ps.data)
    warning("Binomial-response GLMs fitted using the bias-reduction method (brglm)")
  }

  myfit[[1]] = model.glm

  if (ps.link == "probit") {
    myfit[[2]] = PsychDelta(model.glm)
    if (ps.lines == T) {
      segments(x0 = myfit[[2]][1, 3], y0 = 0.5, x1 = myfit[[2]][1, 4], y1 = 0.5, col = ps.col,
               lty = ps.lty, lwd = ps.lwd)
    }
  } else {
    myfit[[2]] = NA
    warning("Use the probit link function to get the estimate of the PSE and the JND")
  }

  if (ps.lines == T) {
    if (is.na(ps.x)) {
      ps.x <- data.frame(seq(x.range[1], x.range[2], length.out = 1000))
      names(ps.x) = attr(ps.terms, "term.labels")
    }
    lines(x = seq(x.range[1], x.range[2], length.out = 1000),
          y = predict(object = model.glm, newdata = data.frame(ps.x), type = "response"),
          col = ps.col, lty = ps.lty, lwd = ps.lwd)
    
  }
  return(myfit)
}

#' Plotting Psychometric Functions given PSE and JND
#'
#' \code{PsychShape()} plots a psychometric function with known pse and jnd
#' on an existing plot.
#'
#' @param pse,jnd the pse and the jnd of the desired psychometric function
#' @param ps.link a link function for the binomial family of error distribution (see Details).
#' @param x.range a vector of length two specifying the range of the function
#' @param ps.col color of the line to be plotted
#' @param ps.lwd line width
#' @param ps.lty line type
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
#' y = c(0,1)
#' x = c(-40, 40)
#' plot(y ~ x, type = "n", bty = "n", lab = c(5,3,7))
#' PsychShape(pse = 0, jnd = 6, x.range = c(-40, 40), ps.col = "gray", ps.lwd = 3)
#' PsychShape(pse = 6, jnd = 6, x.range = c(-40, 40), ps.col = "black")
#' PsychShape(pse = 6, jnd = 6, x.range = c(-40, 40), ps.col = "red", ps.link = "logit", ps.lwd = 3)
#'
#' @importFrom stats plogis
#' @export
#' 
PsychShape <- function(pse = 0, jnd, x.range = c(NA, NA), ps.link = "probit", ps.col = "black", ps.lwd = 1,
                       ps.lty = "solid") {
  x = NA
  if (ps.link == "probit") {
    slope = qnorm(0.75) * (1/jnd)
    curve(expr = pnorm(x, mean = pse, sd = 1/slope), from = x.range[1], to = x.range[2], col = ps.col,
          add = TRUE, lwd = ps.lwd, lty = ps.lty)
  } else {
    if (ps.link == "logit") {
      slope = log(3) * (1/jnd)
      curve(expr = plogis(x, location = pse, scale = 1/slope), from = x.range[1], to = x.range[2],
            col = ps.col, add = TRUE, lwd = ps.lwd, lty = ps.lty)
    } else {
      warning("The function only works with probit and logit link function")
    }
  }
  return(0)
}
