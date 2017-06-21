#' PSE and JND  Using Delta Method
#'
#' Estimates the Point of Subjective Equivalence (PSE), the Just Noticeable
#' Difference (JND) and the related Standard Errors by means of Delta Method.
#' The method only applies to univariable GLMs (psychometric functions) having
#' a \emph{probit} link function.
#'
#' @details \code{PsychDelta} estimates PSE and JND of a univariable psychometric
#' function (object of class \code{"glm"}).
#'
#' @param model the fitted psychometric function. An object of class \code{"glm"}.
#' @param alpha significance level of the confidence interval.
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
#' Faraggi D, Izikson P, Reiser B (2003). Confidence intervals for the
#' 50 per cent response dose. Statistics in Medicine, 22(12), 1977b-1988.
#'
#' Moscatelli A, Mezzetti M, Lacquaniti F (2012). Modelling Psychophysical Data
#' at the Population-Level: The Generalized Linear Mixed Model.
#' Journal of Vision, 12(11):26, 1-17.
#'
#' @seealso \code{\link[stats]{glm}} for for Generalized Linear Models (without
#' random effects) and \code{\link[lme4]{glmr}} for Generalized Linear Mixed
#' Models (including random effects). \code{MERdelta.probit} and \code{MERtreatment}
#' for univarible and multivariable GLMM, respectively (object of class
#' \code{"\linkS4class{merMod}"}). \code{\link{pseMer}} provides the bootstrap-based
#' confidence intervals.
#'
#' @examples
#' #load simulated data
#' data(psych)
#' #fit a glm (probit link)
#' model.glm = glm(formula = cbind(Longer, Total - Longer) ~ X,
#' family = binomial(link = "probit"), data = psych)
#' PsychDelta(model.glm)
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
#' \code{PsychFunction} is used to fit psychometric functions using either
#' \code{glm()} or \code{brglm()}, estimate the PSE, the JND and the related
#' confidence intervals and draw the curve on an existing plot.
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
#' @param ps.lty line type of the lines to be plotted
#' @param br  logical. If TRUE, brglm is used if fitted values are equal to 0 or 1
#'
#' @details If lines is TRUE, the function draw model predictions on an existing plot.
#' Only for univariable glm of the type \code{F(Y) ~ X}, where X is a continuous
#' predictor. If ps.x is empty, the new data frame is a vector of length = 1000,
#' whose range is specified from x.range. Std. Errors and 95\% confidence intervals
#' of the PSE and JND are estimated via Delta Methods, see Faraggi et al (2003).
#'
#' @return \code{PsychFunction} returns list including the fitted glm (or brglm),
#' the estimate of PSE and JND and a flag to indicate if brglm was called.
#'
#' @references
#' Faraggi D, Izikson P, Reiser B (2003). Confidence intervals for the 50 per
#' cent response dose. Statistics in Medicine, 22(12), 1977b-1988.
#'
#' Moscatelli A, Mezzetti M, Lacquaniti F (2012). Modelling Psychophysical Data
#' at the Population-Level: The Generalized Linear Mixed Model.
#' Journal of Vision, 12(11):26, 1-17.
#'
#' @seealso \code{\link[stats]{glm}} for for Generalized Linear Models.
#'
#'
#'# simulate data from a single participant
#'datafr.S1 <- PsySimulate(fixeff = c(-7.5, 0.0875), nsubject = 1, constant = T) 

#' #fit a glm (probit link)
#' model.glm = glm(formula = cbind(Longer, Total - Longer) ~ X,
#' family = binomial(link = "probit"), data = datafr.S1)
#'
#' #psychometric function to fit single-subject data
#' plot(Longer/Total ~ X, data = datafr.S1)
#' fit.S1 = PsychFunction(ps.formula = cbind(Longer, Total - Longer) ~ X,
#'                         ps.link = "probit", ps.data = datafr.S1,
#'                         x.range = c(40, 120), ps.lines = T)
#' @export
#' @importFrom brglm brglm
#' 
PsychFunction <- function(ps.formula, ps.link, ps.data, x.range = c(NA, NA), ps.x = NA, ps.lines = F,
                          ps.col = "black", ps.lty = "dashed", br = F) {
  myfit = vector("list", 3)
  ps.terms = terms(ps.formula)

  model.glm = glm(formula = ps.formula, family = binomial(link = ps.link), data = ps.data)
  brflag = ifelse(trunc(max(model.glm$fitted.values)) == 1 & trunc(min(model.glm$fitted.values)) ==
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
               lty = ps.lty)
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
    lines(x = seq(x.range[1], x.range[2], length.out = 1000), y = predict(object = model.glm,
                                                                          newdata = data.frame(ps.x), type = "response"), col = ps.col, lty = ps.lty)
  }
  return(myfit)
}

#' Fitting and Plotting Psychometric Functions
#'
#' \code{PsychShape()} plots a psychometric function of a known pse and jnd
#' on an existing plot.
#'
#' @param pse,jnd the pse and the jnd of the desired psychometric function
#' @param ps.link a link function for the binomial family of error distribution.
#' See `Details'
#' @param x.range a vector of length two specifying the range of the function
#' @param ps.col color of the line to be plotted
#' @param ps.lwd line width
#' @param ps.lty line type
#'
#' @details PsychShape() can be used to visualize the predicted results of a
#' psychophysical experiment or to plot a fitted psychometric function whose
#' values of pse and jnd are known. Currently only working with probit and logit
#' link function.
#'
#' @references
#' Moscatelli A, Mezzetti M, Lacquaniti F (2012). Modelling Psychophysical Data
#' at the Population-Level: The Generalized Linear Mixed Model.
#' Journal of Vision, 12(11):26, 1-17.
#'
#' @seealso \code{\link[stats]{glm}} for for Generalized Linear Models.
#'
#' @examples
#' y = c(0,1)
#' x = c(-40, 40)
#' plot(y ~ x, type = "n", bty = "n", lab = c(5,3,7))
#' PsychShape(pse = 0, jnd = 6, x.range = c(-40, 40), ps.col = "gray", ps.lwd = 3)
#' PsychShape(pse = 6, jnd = 6, x.range = c(-40, 40), ps.col = "black")
#' PsychShape(pse = 6, jnd = 6, x.range = c(-40, 40), ps.col = "red", ps.link = "logit", ps.lwd = 3)
#'
#' @export
PsychShape <- function(pse = 0, jnd, x.range = c(NA, NA), ps.link = "probit", ps.col = "black", ps.lwd = 1,
                       ps.lty = "solid") {
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
