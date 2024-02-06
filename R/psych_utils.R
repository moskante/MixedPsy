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
#' Journal of Vision, 12(11):26, 1-17. doi:10.1167/12.11.26
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
#' @importFrom stats vcov qnorm
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
#' Journal of Vision, 12(11):26, 1-17. doi:10.1167/12.11.26
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
#' psych.S1 <- PsychFunction(formula = cbind(Longer, Total - Longer) ~ X, 
#' link = "probit", data = data.S1)
#' plotP1 <- PsychPlot(psych.S1$glm, showData = TRUE, ps.lab = "S1") 
#' 
#' data.S2 <- subset(simul_data, Subject == "S2")
#' glm.S2 <- glm(formula = cbind(Longer, Total - Longer) ~ X, 
#'            family = binomial(link = "probit"), data = data.S2)
#' plotP2 <- PsychPlot(glm.S2, addTo = plotP1, ps.lab = "S2")
#'
#' @keywords GLM plot
#' 
#' @import ggplot2
#' @importFrom stats predict
#' 
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
#' Journal of Vision, 12(11):26, 1-17. doi:10.1167/12.11.26
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
#' @importFrom stats qnorm pnorm plogis
#' 
#' @export
#' 
#' @keywords GLM plot
#' 
PsychShape <- function(pse = 0, jnd = 1, p = 0.75, x.range = c(NA, NA), ps.link = c("probit"), 
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


#' Fit Psychometric Function from Single-Subject Response
#'
#' This function provides an interface for fitting psychometric functions using various models: generalized linear model (glm), glm with bias reduction (brglm), generalized non-linear model (gnlm).
#'
#' @param formula A formula specifying the model (glm or brglm). If NULL, 'response' and 'stimuli' must be provided. 
#' @param response A character vector with the name of the response variables. It is used to create the matrix of the binomial response. It must be provided for gnlm.
#' @param stimuli A character vector with the name of the stimulus variable. It must be provided for gnlm.
#' @param model A character string specifying the model. Possible options are 'glm' for the generalized linear model, 'brglm' for glm with bias reduction, 'gnlm' for generalized non-linear model. Default is 'glm'.
#' @param link A character string specifying the link function. For generalized linear models, it defines the link function specified in the family object. See \code{\link[stats]{family}} for details of family functions (for binomial family). For generalized non-linear models, it defines the function being fitted. Possible options are 'probit' for the cumulative normal distribution, 'logit' for the cumulative logit distribution, 'weibull' for cumulative Weibull distribution.
#' @param data A data frame containing the variables used in the model.
#' @param guess Logical or numeric value indicating whether to include a guessing parameter. Only used if model = 'gnlm'. Default is FALSE.
#' @param lapse Logical or numeric value indicating whether to include a lapse parameter. Only used if model = 'gnlm'. Default is FALSE.
#'
#' @return A list containing the fitted models and additional information. 
#'
#' @examples
#' \dontrun{
#' PsychFunction_new(formula = y ~ x, data = my_data)
#' }
#' @export
PsychFunction <- function (formula = NULL, response = NULL, stimuli = NULL, model = "glm", link = "probit", data, guess = FALSE, lapse = FALSE){
  
  stopifnot(is.character(model), is.character(link), is.data.frame(data))
  
  myfit = list()
  
  model_glm <- PsychFunction_glm(formula, response, stimuli, model = "glm", link, data)
  
  myfit$glm <- model_glm$model
  myfit$recommend_br <- model_glm$flag    
  
  if(model == "brglm"){
    model_brglm <- PsychFunction_glm(formula, response, stimuli, model, link, data)
    myfit$brglm <- model_brglm$model
  }else if(model == "gnlm"){
    x_values <- data[,stimuli]
    myfit$gnlm <- PsychFunction_gnlm(myfit$glm, link, response, data, stimuli, guess, lapse)
  }
  
  return(myfit)
}




#' Internal Function: Fit Generalized Linear Models
#'
#' This function fits a generalized linear model for psychometric functions using glm or brglm.
#'
#' @param formula A formula specifying the model. If NULL, 'response' and 'stimuli' must be provided.
#' @param response A character vector of response variables.
#' @param stimuli A character vector of stimulus variables.
#' @param model A character string specifying the model ('glm', 'brglm').
#' @param link A character string specifying the link function ('probit', 'logit', 'weibull').
#' @param data A data frame containing the variables specified in the formula.
#'
#' @return A list containing the fitted model and additional information.
#'
#' @importFrom brglm brglm
#' @importFrom stats glm
#' @importFrom stats as.formula
#' 
PsychFunction_glm <- function(formula, response, stimuli, model, link, data){
  if (is.null(formula)){
    stopifnot(is.character(response), is.character(stimuli))
    formula_string <- paste("cbind(", response[1], ",", response[2], ")", "~", stimuli)
    formula <- as.formula(formula_string)  
  }
  if (link == "weibull"){
    print("Warning: weibull is not a possible link function for glm. probit was used instead")
    link = "probit"
  }
  model_glm <- glm(formula, family = binomial(link = link), 
                   data = data)
  
  eps <- 1e-15
  brflag <- ifelse(1 - max(model_glm$fitted.values) <= eps & 
                     trunc(min(model_glm$fitted.values)) == 0, TRUE, FALSE)
  if (model == "brglm" ) { #& brflag == TRUE
    model_glm <- brglm(formula, family = binomial(link = link), 
                       data = data)
  }
  
  return(list(model = model_glm, flag = brflag))
}

#' Internal Function: Fit Generalized Nonlinear Models
#'
#' This function fits a generalized nonlinear model for psychometric functions.
#'
#' @param model_glm A glm object obtained from PsychFunction_glm.
#' @param link A character string specifying the link function ('probit', 'logit', 'weibull').
#' @param response A character vector of response variables.
#' @param stimuli A character vector of stimulus variables.
#' @param data A data frame containing the variables specified in the formula.
#' @param guess Logical or numeric value indicating whether to include a guessing parameter.
#' @param lapse Logical or numeric value indicating whether to include a lapse parameter.
#' @importFrom gnlm gnlr
#' @importFrom here here
#' @importFrom rmutil fnenvir
#' @return A gnlr object representing the fitted model.
#'
PsychFunction_gnlm <- function(model_glm, link, response, data, stimuli, guess, lapse){
  #source(here("R", "global.R"))
  # to do: warning if formula is defined instead of response and stimuli - or get response and stimuli from formula but give a warning.
  response <- parse(text = paste("cbind(", response[1], ",", response[2], ")"))
  
  
  glm_coeff <- summary(model_glm)$coefficients[,1]
  start_estimate <- c(-glm_coeff[1]/glm_coeff[2], ifelse(link == "weibull", 2, 1/glm_coeff[2]))
  
  lambda <- if (isTRUE(lapse)) runif(1, min = 0, max = 0.05) else if (is.numeric(lapse)) lapse else FALSE
  gamma <- if (isTRUE(guess)) runif(1, min = 0, max = 0.05) else if (is.numeric(guess)) guess else FALSE
  
  if (is.numeric(gamma) && !is.numeric(lambda)){
    start_estimate <- c(start_estimate, gamma)
  }else if (!is.numeric(gamma) && is.numeric(lambda)){
    start_estimate <- c(start_estimate, lambda)
  }else if(is.numeric(gamma) && is.numeric(lambda)){
    start_estimate <- c(start_estimate, gamma, lambda)
  }
  
  #x_values <- data[,stimuli]
  #this_env <- new.env()
  #x_values <- data[,stimuli]
  #environment(x_values) <- NULL
  assign_global_var(data,stimuli)
  switch_function <- fnenvir(function(p) pnorm(x_values, mean = p[1], sd = p[2]))
  #switch_mu_function(x_values, func_name = link, data, stimuli, gamma, lambda)
  #attr(switch_function, "model") <- TRUE
  #environment(switch_function) <- this_env
  model_gnlm <- gnlr(y = with(data, eval(response)), distribution = "binomial",
                     mu = switch_function, pmu = start_estimate)
  
  #rmGlobalVar()
  return(model_gnlm)
}

assign_global_var <- function(data, stimuli) {
  x_values <<- data[,stimuli]
}


#' Internal Function: Switch mu Function
#'
#' This function switches between different mu functions based on the provided parameters.
#'
#' @param func_name A character string specifying the mu function ('probit', 'logit', 'weibull').
#' @param data dataset
#' @param stimuli name of stimulus variable
#' @param gamma A numeric or logical value indicating the gamma parameter.
#' @param lambda A numeric or logical value indicating the lambda parameter.
#'
#' @return A function representing the selected mu function.
#' @importFrom stats pweibull
#' 
switch_mu_function <- function(x_values, func_name, data, stimuli, gamma, lambda) {
  #setGlobalVar(data, stimuli) #x_values defined as global variable (<<-) due to gnlm syntax. 
  #x_values <- parent.frame()$x_values
  function(p) pnorm(x_values, mean = p[1], sd = p[2])
  
  # if (isFALSE(gamma) && isFALSE(lambda)){
  #   #switch(func_name,
  #   #       probit = function(p) pnorm(x_values, mean = p[1], sd = p[2]),
  #   #       logit = function(p) plogis(x_values, location = p[1], scale = p[2]),
  #   #       weibull = function(p) pweibull(x_values, scale = p[1], shape = p[2]))
  # }else if (is.numeric(gamma) && isFALSE(lambda)){
  #   switch(func_name,
  #          probit = function(p) p[3] + (1 - p[3]) * pnorm(x_values, mean = p[1], sd = p[2]),
  #          logit = function(p) p[3] + (1 - p[3]) * plogis(x_values, location = p[1], scale = p[2]),
  #          weibull = function(p) p[3] + (1 - p[3]) * pweibull(x_values, scale = p[1], shape = p[2]))
  # }else if (isFALSE(gamma) && is.numeric(lambda)){
  #   switch(func_name,
  #          probit = function(p) (1 - p[3]) * pnorm(x_values, mean = p[1], sd = p[2]),
  #          logit = function(p) (1 - p[3]) * plogis(x_values, location = p[1], scale = p[2]),
  #          weibull = function(p) (1 - p[3]) * pweibull(x_values, scale = p[1], shape = p[2]))
  # }else if(is.numeric(gamma) && is.numeric(lambda)){
  #   switch(func_name,
  #          probit = function(p) p[3] + (1 - p[3] - p[4]) * pnorm(x_values, mean = p[1], sd = p[2]),
  #          logit = function(p) p[3] + (1 - p[3] - p[4]) * plogis(x_values, location = p[1], scale = p[2]),
  #          weibull = function(p) p[3] + (1 - p[3] - p[4]) * pweibull(x_values, scale = p[1], shape = p[2]))
  # }
}