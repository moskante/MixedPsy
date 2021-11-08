#' PSE/JND from GLM Using Delta Method
#'
#' Estimate Point of Subjective Equivalence (PSE), Just Noticeable
#' Difference (JND), and related Standard Errors of an individual participant 
#' by means of Delta Method.
#' The method only applies to a GLM (object of class \code{\link[stats]{glm}}) with one continuous 
#' predictor and a \emph{probit} link function.
#'
#' @param model.obj the fitted psychometric function. An object of class \code{\link[stats]{glm}}.
#' @param alpha significance level of the confidence interval.Default is 0.05 (95\% confidence interval).
#' @param p probability value relative to the JND upper limit. Default is 0.75 (value for 50\% JND).
#' 
#' @details \code{PsychDelta} estimates PSE and JND of a psychometric
#' function (object of class \code{glm}).
#' 
#' @return \code{PsychDelta} returns a matrix including estimate, standard error,
#' inferior and superior bounds of the confidence interval of PSE and JND. Confidence Intervals
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
#' Knoblauch, K., & Maloney, L. T. (2012). Modeling psychophysical data in R (Vol. 32). 
#' Springer Science & Business Media.
#'
#' Moscatelli, A., Mezzetti, M., & Lacquaniti, F. (2012). Modeling psychophysical data 
#' at the population-level: The generalized linear mixed model. 
#' Journal of Vision, 12(11):26, 1-17. https://doi.org/10.1167/12.11.26
#'
#' @seealso \code{\link[stats]{glm}} for fitting a Generalized Linear Model to a single-subject response. \code{\link[lme4]{glmer}} 
#' for Generalized Linear Mixed Models (including fixed and random effects). \code{MixDelta} for estimating PSE and JND at a population level 
#' with delta method. 
#' 
#' @keywords DeltaMethod GLM 
#'
#' @examples
#' data.S1 <- subset(simul_data, Subject == "S1")
#' model.glm = glm(formula = cbind(Longer, Total - Longer) ~ X,
#' family = binomial(link = "probit"), data = data.S1)
#' PsychDelta(model.glm)
#' 
#' @export
#'
PsychDelta <- function(model.obj, alpha = 0.05, p = 0.75) {

    pse <- -model.obj$coef[1]/model.obj$coef[2]
    BETA <- model.obj$coef[2]

    var.alpha <- vcov(model.obj)[1, 1]
    var.beta <- vcov(model.obj)[2, 2]
    cov.alpha.beta <- vcov(model.obj)[2, 1]

    var.pse <- (1/BETA^2) * (var.alpha + (2 * pse * cov.alpha.beta) + (pse^2 * var.beta))  #PSE
    inferior.pse <- pse - (qnorm(1 - (alpha/2)) * sqrt(var.pse))
    superior.pse <- pse + (qnorm(1 - (alpha/2)) * sqrt(var.pse))

    jnd <- qnorm(p) * (1/BETA)
    var.jnd <- (qnorm(p) * (-1/BETA^2))^2 * var.beta  #JND
    inferior.jnd <- jnd - (qnorm(1 - (alpha/2)) * sqrt(var.jnd))
    superior.jnd <- jnd + (qnorm(1 - (alpha/2)) * sqrt(var.jnd))

    output <- matrix(rbind(c(pse, sqrt(var.pse), inferior.pse, superior.pse), c(jnd, sqrt(var.jnd),
        inferior.jnd, superior.jnd)), nrow = 2, dimnames = list(param <- c("pse", "jnd"), statistics <- c("Estimate",
        "Std. Error", "Inferior", "Superior")))

    return(output)
}


#' Psychometric Function and PSE/JND Parameters from Single-Subject Response
#'
#' Fit psychometric functions using \code{\link[stats]{glm}} or \code{\link[brglm]{brglm}}. 
#' Estimate PSE, JND, and related confidence intervals with Delta Method. 
#'
#' @param ps.formula an object of class \code{\link[stats]{formula}}, such as \code{cbind(yes, no) ~ X}
#' @param ps.link link function for the binomial family of error distribution. Default is \code{probit}.
#' @param ps.data a data frame including the variables used in the model.
#' @param br  logical. If TRUE, \code{\link[brglm]{brglm}} for bias reduction is used if values are equal to 0 or 1. 
#' Default is FALSE.
#'
#' @details Estimates are computed only for GLM of the type \code{F(Y) ~ X}, where X is a continuous
#' predictor. Std. Errors and 95\% confidence intervals of PSE and JND are estimated via Delta Methods. 
#' Currently only working with \emph{probit} link function.
#'
#' @return \code{PsychFunction} returns a list including the fitted model,
#' the estimate of PSE and JND and a flag to indicate if \code{\link[brglm]{brglm}} was called.
#'
#' @note \code{PsychFunction} returns the same parameter estimate as \code{\link{PsychDelta}}, without an explicit call to \code{\link[stat]{glm}}. 
#' Moreover, it allows to fit the model using \code{\link[brglm]{brglm}} in case of complete or quasi separation.
#' 
#' @references
#' Faraggi, D., Izikson, P., & Reiser, B. (2003). Confidence intervals for the 50 per cent 
#' response dose. Statistics in medicine, 22(12), 1977-1988. https://doi.org/10.1002/sim.1368
#'
#' Moscatelli, A., Mezzetti, M., & Lacquaniti, F. (2012). Modeling psychophysical data 
#' at the population-level: The generalized linear mixed model. 
#' Journal of Vision, 12(11):26, 1-17. https://doi.org/10.1167/12.11.26
#' 
#' @seealso \code{\link[stats]{glm}} for Generalized Linear Models. 
#' \code{\link[brglm]{brglm}} for fitting a GLM using bias reduction.
#' \code{\link{PsychPlot}} for plotting a psychometric function given a \code{\link[stats]{glm}} (or \code{\link[brglm]{brglm}}) object.
#' \code{\link{PsychPlot}} for plotting a a psychometric function from a GLM. 
#' \code{\link{PsychShape}} for plotting a psychometric function given PSE and JND. 
#' 
#' @keywords GLM DeltaMethod
#'
#' @examples
#' data.S1 <- subset(simul_data, Subject == "S1")
#' psych.S1 <- PsychFunction(ps.formula = cbind(Longer, Total - Longer) ~ X, 
#' ps.link = "probit", ps.data = data.S1)
#'                         
#' @importFrom brglm brglm
#' @export
#'
<<<<<<< HEAD
PsychFunction <-  function (ps.formula, ps.link, ps.data, br = F) {
  
  myfit = list()
=======
PsychFunction <- function(ps.formula, ps.link, ps.data, x.range = c(NA, NA), ps.x = NA, ps.lines = F,
                          ps.col = "black", ps.lty = "dashed", ps.lwd = 1, br = F) {
  myfit = vector("list", 3)
>>>>>>> master
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
<<<<<<< HEAD
    myfit$estimate = PsychDelta(model.glm)
=======
    myfit[[2]] = PsychDelta(model.glm)
    if (ps.lines == T) {
      segments(x0 = myfit[[2]][1, 3], y0 = 0.5, x1 = myfit[[2]][1, 4], y1 = 0.5, col = ps.col,
               lty = ps.lty, lwd = ps.lwd)
    }
>>>>>>> master
  } else {
    myfit$estimate = NA
    warning("Use the probit link function to get the estimate of the PSE and the JND")
  }
  
  myfit$info <- list(brflag = brflag)
  
  return(myfit)
}

<<<<<<< HEAD
#' Plot Psychometric Function from GLM
#'
#' Plot a psychometric function given an object of class \code{\link[stats]{glm}} or \code{\link[brglm]{brglm}}. 
#' The plot can be drawn on a new or existing \code{ggplot} object.
#'
#' @param model.obj the fitted psychometric function. An object of class \code{\link[stats]{glm}} or \code{\link[brglm]{brglm}}.
#' @param addTo specifies an existing \code{ggplot} object where the new line should be plotted. If no object is given, the function is drawn on a new plot.
#' @param showData logical, defines if proportion of binomial responses for each stimulus level are presented. Default is TRUE.
#' @param ps.type,ps.size type and size of the plotted line (see \code{"ggplot2-spec"}).
#' @param ps.lab label assigned to the psychometric curve. The label is coded by the color aesthetic. 
#'
#' @return \code{PsychPlot} returns a \code{\link[ggplot2]{ggplot}} object. 
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
#' \code{\link{MixPlot}} for plotting individual responses from a GLMM.
#'
#' @examples
#' data.S1 <- subset(simul_data, Subject == "S1")
#' psych.S1 <- PsychFunction(ps.formula = cbind(Longer, Total - Longer) ~ X, 
#' ps.link = "probit", ps.data = data.S1)
#' plotP1 <- PsychPlot(psych.S1$model, showData = TRUE, ps.lab = "S1") 
#' 
#' data.S2 <- subset(simul_data, Subject == "S2")
#' glm.S2 <- glm(formula = cbind(Longer, Total - Longer) ~ X, 
#'            family = binomial(link = "probit"), data = data.S2)
#' plotP2 <- PsychPlot(glm.S2, addTo = plotP1, ps.lab = "S2")
#'
#' @keywords GLM plot
#' 
#' @import ggplot2
#' @export
#' 
PsychPlot <- function(model.obj, addTo = NULL, showData = TRUE,
                      ps.type = "solid", ps.size = 1, ps.lab = ""){
  if(is.null(addTo)){
    p <- ggplot()
  }else{ 
    p <- addTo}
  
  # GLM
  temp <- names(model.obj$model)
  xname <- temp[2]
  
  data <- model.obj$data
  data$y <- model.obj$y
  longData <- data.frame(pretty(data[[xname]], 1000))
  names(longData) <- xname
  longData$y <- predict(object = model.obj, newdata = longData, type = "response")
  
  plot <- list()
  
  plot$ps <- geom_line(inherit.aes = FALSE, data = longData, aes_string(xname, 'y', color = 'ps.lab'), 
                       linetype = ps.type, size = ps.size)  
  if(isTRUE(showData)){
    plot$data <- geom_point(inherit.aes = FALSE, data = data, aes_string(xname, 'y', color = 'ps.lab'))
  }
  
  if(model.obj$family[["link"]] == "probit"){
    psych <- PsychDelta(model.obj)
    plot$segment <- geom_segment(inherit.aes = FALSE, data = NULL, 
                                 aes(x = psych["pse", "Inferior"], xend = psych["pse", "Superior"], 
                                     y = 0.5, yend = 0.5, color = ps.lab),
                                 size = ps.size)
  } 
  print(p + plot)
  
=======
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
>>>>>>> master
}

#' Plot Psychometric Functions given PSE and JND
#'
#' Plot a psychometric function with known PSE and JND on a new or existing \code{ggplot} object.
#'
#' @param pse,jnd point of subjective equivalende (PSE) and just noticeable difference (JND) of the desired psychometric function.
#' @param p probability value relative to the JND upper limit. Default is 0.75 (value for 50\% JND).
#' @param x.range vector of length two specifying the range of the psychometric function.
#' @param ps.link a link function for the binomial family of error distribution.
#' @param ps.type,ps.size,ps.color type, size, and color of the plotted line (see \code{"ggplot2-spec"}).
#' @param addTo specifies an existing \code{ggplot} object where the new line should be plotted. If no object is given, the function is drawn on a new plot.
#'
#' @details \code{PsychShape()} can be used to visualize the predicted results of a
#' psychophysical experiment or to plot a fitted psychometric function whose
#' values of pse and jnd are known. Currently only working with probit and logit
#' link function.
#' 
#' @return \code{PsychShape} returns a \code{\link[ggplot2]{ggplot}} object. 
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
#' \code{\link{PsychFunction}} and \code{\link{PsychDelta}} for estimation of PSE and JND from response data.
#' \code{\link{PsychPlot}} for plotting a a psychometric function from a GLM. 
#'
#' @examples
#' p <- PsychShape(pse = 0, jnd = 6, x.range = c(-40, 40), ps.color = "gray", ps.size = 3)
#' p1 <- PsychShape(pse = 6, jnd = 6, x.range = c(-40, 40), ps.col = "black", addTo = p)
#' p2 <- PsychShape(pse = 6, jnd = 6, x.range = c(-40, 40), ps.col = "red", ps.link = "logit", 
#' ps.type = "dashed", addTo = NULL)
#'
#' @import ggplot2
#' @export
#' 
#' @keywords GLM plot
#' 
PsychShape <- function(pse = 0, jnd = 1, p = 0.75, x.range = c(NA, NA), ps.link = c("probit", "logit"), 
                       ps.type = "solid", ps.size = 1, ps.color = "black", addTo = NULL) {
  if(is.null(addTo)){
    pl <- ggplot()
  }else{ 
    pl <- addTo}
  
  x = pretty(x.range, 100)
  if (ps.link == "probit") {
    slope = qnorm(p) * (1/jnd)
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
  
  print(pl + plot)
}

