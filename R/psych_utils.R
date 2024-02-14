#' PsychEstimate - Estimate Psychological Statistics
#'
#' This function returns Point of Subjective Equivalence (PSE), Just Noticeable
#' Difference (JND), and related Standard Errors of an individual participant using 
#' either the delta method or bootstrapping.
#'
#' @param model_obj An object of class 'gnlm' or 'glm' representing the fitted model.
#' @param method The estimation method. Options are "delta" or "boot".
#' @param ... description
#' @inheritParams PsychDelta
#' @inheritParams PsychBoot
#' 
#' @return This function returns estimates of psychometric parameters.
#' 
#' @examples
#' data.S1 <- subset(simul_data, Subject == "S1")
#' model.glm = glm(formula = cbind(Longer, Total - Longer) ~ X,
#' family = binomial(link = "probit"), data = data.S1)
#' PsychEstimate(model.glm, method = "delta")
#' 
#' @export
#' 
PsychEstimate <- function(model_obj, alpha = 0.05, p = 0.75, method, ...){
  allowed_methods <- c("delta", "boot")
  method <- match.arg(method, allowed_methods)
  
  model_class <- class(model_obj)
  
  if (any(model_class == "gnlm")) {
    if (method == "boot") {
      PsychBoot(model_obj, alpha, p)
    } else {
      stop("Delta method not available for model of class gnlm.")
    }
  } else if (any(model_class == "glm")) {
    if (method == "boot") {
      PsychBoot()
    } else {
      PsychDelta(model_obj, alpha, p)
    }
  } else {
    stop("Unsupported model class.")
  }
  
}

#' Internal function: PSE/JND from GLM Using Delta Method
#'
#' Estimate Point of Subjective Equivalence (PSE), Just Noticeable
#' Difference (JND), and related Standard Errors of an individual participant 
#' by means of Delta Method.
#' The method only applies to a GLM (object of class \code{\link[stats]{glm}}) with one continuous 
#' predictor and a \emph{probit} or \emph{logit} link function.
#'
#' @param model_obj the fitted psychometric function. An object of class \code{\link[stats]{glm}}.
#' @param alpha significance level of the confidence interval.Default is 0.05 (95\% confidence interval).
#' @param p probability value relative to the JND upper limit. Default is 0.75 (value for 50\% JND).
#' 
#' @details \code{PsychDelta} estimates PSE and JND of a psychometric
#' function for an object of class \code{glm}.
#' 
#' @return \code{PsychDelta} returns a matrix including estimate, standard error,
#' inferior and superior bounds of the confidence interval of PSE and JND. Confidence Intervals
#' are computed as: \eqn{Estimate +/- z(1-(\alpha/2)) * Std.Error}.
#'
#' @note The function assumes that the first model coefficient is the intercept
#' and the second is the slope. The estimate of the JND assumes a \emph{probit}
#' or \emph{logit} link function.
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
#' 
#' @importFrom stats vcov qnorm
#' 
#' 
PsychDelta <- function(model_obj, alpha = 0.05, p = 0.75) {
  
  pse <- -model_obj$coef[1]/model_obj$coef[2]
  BETA <- model_obj$coef[2]
  
  var_alpha <- vcov(model_obj)[1, 1]
  var_beta <- vcov(model_obj)[2, 2]
  cov_alpha_beta <- vcov(model_obj)[2, 1]
  
  var_pse <- (1/BETA^2) * (var_alpha + (2 * pse * cov_alpha_beta) + (pse^2 * var_beta))  #PSE
  inferior_pse <- pse - (qnorm(1 - (alpha/2)) * sqrt(var_pse))
  superior_pse <- pse + (qnorm(1 - (alpha/2)) * sqrt(var_pse))
  
  family <- model_obj$family
  if (family$link == "probit"){
    jnd_const <- qnorm(p)
  }else if(family$link == "logit"){
    jnd_const <- log(3)
  }else{
    warning("link not found")
  }
  jnd <- jnd_const * (1/BETA)
  var_jnd <- (jnd_const * (-1/BETA^2))^2 * var_beta  #JND
  inferior_jnd <- jnd - (qnorm(1 - (alpha/2)) * sqrt(var_jnd))
  superior_jnd <- jnd + (qnorm(1 - (alpha/2)) * sqrt(var_jnd))
  
  output <- matrix(rbind(c(pse, sqrt(var_pse), inferior_pse, superior_pse), 
                         c(jnd, sqrt(var_jnd), inferior_jnd, superior_jnd)), 
                   nrow = 2, 
                   dimnames = list(param <- c("pse", "jnd"), 
                                   statistics <- c("Estimate", "Std. Error", "Inferior", "Superior")))
  
  return(output)
}

#' Internal function: PSE/JND from GLM or GNLM with Bootstrapping 
PsychBoot <- function(){
  
}


#' Plot Psychometric Function from GLM
#'
#' Plot a psychometric function given an object of class \code{\link[stats]{glm}} or \code{\link[brglm]{brglm}}. 
#' The plot can be drawn on a new or existing \code{ggplot} object.
#'
#' @param model_obj the fitted psychometric function. An object of class \code{\link[stats]{glm}} or \code{\link[brglm]{brglm}}.
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
PsychPlot <- function(model_obj, addTo = NULL, showData = TRUE,
                      ps.type = "solid", ps.size = 1, ps.lab = ""){
  if(is.null(addTo)){
    p <- ggplot()
  }else{ 
    p <- addTo}
  
  # GLM
  temp <- names(model_obj$model)
  xname <- temp[2]
  
  data <- model_obj$data
  data$y <- model_obj$y
  longData <- data.frame(pretty(data[[xname]], 1000))
  names(longData) <- xname
  longData$y <- predict(object = model_obj, newdata = longData, type = "response")
  
  plot <- list()
  
  plot$ps <- geom_line(inherit.aes = FALSE, data = longData, aes_string(xname, 'y', color = 'ps.lab'), 
                       linetype = ps.type, size = ps.size)  
  if(isTRUE(showData)){
    plot$data <- geom_point(inherit.aes = FALSE, data = data, aes_string(xname, 'y', color = 'ps.lab'))
  }
  
  if(model_obj$family[["link"]] %in% c("probit", "logit")){
    psych <- PsychDelta(model_obj)
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
#' @param formula A formula specifying the linear model (glm or brglm). If NULL, 'response' and 'stimuli' must be provided. 
#' @param response A character vector with the name of the response variables. It is used to create the matrix of the binomial response. It must be provided for gnlm.
#' @param stimuli A character string with the name of the stimulus variable. It must be provided for gnlm.
#' @param model A character string specifying the model. Possible options are 'glm' for the generalized linear model, 'brglm' for glm with bias reduction, 'gnlm' for generalized non-linear model. Default is 'glm'.
#' @param link A character string specifying the link function. For generalized linear models, it defines the link function specified in the family object. See \code{\link[stats]{family}} for details of family functions for the binomial family. For generalized non-linear models, it defines the function being fitted. Possible options are 'probit' for the cumulative normal distribution, 'logit' for the cumulative logit distribution, 'weibull' for cumulative Weibull distribution.
#' @param data A data frame containing the variables used in the model.
#' @param guess,lapse Logical or numeric values indicating whether to include guessing and lapse parameters, respectively. These parameters are only used if \code{model = 'gnlm'}. Default is FALSE for both. If parameters are FALSE, they are not included in the model. If TRUE, parameters are estimated with a randomly assigned starting value. If numeric, the value is used as a starting estimate. 
#' @param ... Additional arguments included for back compatibility. Deprecated arguments from older package versions ('ps.formula', 'ps.link', 'ps.data', 'br') will be automatically handled.
#'
#' @return A list containing:
#' \item{glm}{The fitted GLM model.}
#' \item{recommend_br}{A logical value indicating whether bias reduced is recommended.}
#' \item{brglm}{The fitted Bias Reduced GLM, if \code{model = "brglm"}.}
#' \item{gnlm}{The fitted generalized non linear model (gnlm), if \code{model = "gnlm"}.}
#' \item{gnlm_coeff}{The estimated coefficients of the gnlm. When guess and lapse parameters are present, the fitted coefficients returned by the \code{glnr} function undergo a log-transformation. This transformation is applied as a result of fitting the parameters of the model using the exponential function, which ensures that the coefficients remain positive.}
#' 
#' @note
#' If 'model' is specified as "gnlm", a generalized linear model (glm) will be fitted to extract initial estimates for the parameters of the non-linear model. 
#' 
#' 
#' @examples
#' \dontrun{
#' # Fit a psychometric function using glm with default probit link function
#' PsychFunction(formula = cbind(n_yes, n_no) ~ x, data = my_data)
#' 
#' # Fit a psychometric function using brglm with link function
#' PsychFunction(formula = cbind(n_yes, n_no) ~ x, model = "brglm", link = "logit", data = my_data)
#' 
#' # Fit a psychometric function using gnlm with guessing and lapse parameters
#' PsychFunction(response = c("n_yes", "n_no"), stimuli = "x", model = "gnlm", 
#'               link = "weibull", guess = TRUE, lapse = TRUE, data = my_data)
#' }
#'
#' @seealso
#' \code{\link{glm}}, \code{\link{brglm}}, \code{\link{gnlr}}
#'
#' @export
PsychFunction <- function(formula = NULL, response = NULL, stimuli = NULL, model = "glm", link = "probit", data, guess = FALSE, lapse = FALSE, ...){
  
  allowed_models <- c("glm", "brglm", "gnlm")
  model <- match.arg(model, allowed_models)
  
  allowed_links <- c("probit", "logit", "weibull")
  link <- match.arg(link, allowed_links)

  handle_deprecated_argument <- function(old_name, new_name) {
    if (old_name %in% names(list(...))) {
      message(paste("In PsychFunction: the '", old_name, "' argument is deprecated. Please use '", new_name, "' instead.", sep = ""))
      return(list(...)[[old_name]])
    }
  }
  
  old_arguments <- c("ps.formula", "ps.link", "ps.data", "br")
  # Check for additional arguments and raise an error if found
  additional_args <- setdiff(names(list(...)), old_arguments)
  if (length(additional_args) > 0) {
    stop(paste("Additional arguments found:", paste(additional_args, collapse = ", ")))
  }
  
  formula <- if ("ps.formula" %in% names(list(...))) handle_deprecated_argument("ps.formula", "formula") else formula
  link <- if ("ps.link" %in% names(list(...))) handle_deprecated_argument("ps.link", "link") else link
  data <- if ("ps.data" %in% names(list(...))) handle_deprecated_argument("ps.data", "data") else data
  if ("br" %in% names(list(...))) {
    message("In PsychFunction: 'br' was removed. For bias reduction glm, use model = 'brglm'")
  }
  
  stopifnot(is.character(model), is.character(link), is.data.frame(data))
  
  myfit = list()
  
  model_glm <- PsychFunction_glm(formula, response, stimuli, model = "glm", link, data)
  
  myfit$glm <- model_glm$model
  myfit$recommend_br <- model_glm$flag    
  
  if(model == "brglm"){
    model_brglm <- PsychFunction_glm(formula, response, stimuli, model, link, data)
    myfit$brglm <- model_brglm$model
  }else if(model == "gnlm"){
    x_values<- data[,stimuli]
    this_envir = parent.frame()
    model_gnlm <- PsychFunction_gnlm(myfit$glm, response, stimuli, link, data, guess, lapse)
    myfit$gnlm <- model_gnlm$model
    myfit$gnlm_coeff <- model_gnlm$gnlr_coeff
  }
  
  return(myfit)
}




#' Internal Function: Fit Generalized Linear Models
#'
#' This function fits a generalized linear model for psychometric functions with glm() or brglm().
#'
#' @param formula A formula specifying the model. If NULL, 'response' and 'stimuli' must be provided.
#' @param response A character vector of response variables.
#' @param stimuli A character string of the stimulus variable.
#' @param model A character string specifying the model ('glm' or 'brglm').
#' @inheritParams PsychFunction
#'
#' @return A list containing:
#' \item{model}{The fitted model.}
#' \item{flag}{A logical value indicating whether the a Bias Reduced GLM (brglm) is recommended.}
#'
#' @importFrom brglm brglm
#' @importFrom stats glm
#' @importFrom stats as.formula
#' 
#' @seealso
#' [\code{\link{PsychFunction}} ]
#' 
PsychFunction_glm <- function(formula, response, stimuli, model, link, data){
  if(is.null(formula)){
    stopifnot(is.character(response), is.character(stimuli))
    formula_string <- paste("cbind(", response[1], ",", response[2], ")", "~", stimuli)
    formula <- as.formula(formula_string) 
    message(paste("glm formula was built from response and stimuli strings as: ", formula_string))
  }else{
    message("The provided formula was used as an argument in glm()")
  }
  
  if (link == "weibull"){
    warning("Weibull is not a possible link function for glm. Probit link was used for fitting the model.")
    link = "probit"
  }
  
  model_glm <- glm(formula, family = binomial(link = link), 
                   data = data)
  
  eps <- 1e-15
  brflag <- ifelse(1 - max(model_glm$fitted.values) <= eps & 
                     trunc(min(model_glm$fitted.values)) == 0, TRUE, FALSE)
  
  if (model == "brglm" ) {
    model_glm <- brglm(formula, family = binomial(link = link), data = data)
  }
  
  return(list(model = model_glm, flag = brflag))
}

#' Internal Function: Fit Generalized Nonlinear Models
#'
#' This function fits a generalized nonlinear model for psychometric functions.
#'
#' @param model_glm A glm object obtained from PsychFunction_glm. This is used to extract starting estimate of the parameters of the non-linear model.
#' @param response A character vector with the name of the response variables.
#' @param stimuli A string with name of the stimuli variable.
#' @inheritParams PsychFunction
#' 
#' @importFrom gnlm gnlr
#' @importFrom here here
#' 
#' @seealso 
#' [\code{\link{PsychFunction}}] [\code{\link{switch_mu_function}} ]
#' 
#' @return A list containing:
#' \item{model_gnlm}{The fitted generalized non-linear model.}
#' \item{gnlr_coeff}{The estimated coefficients (see notes). }
#' 
#' @note For easier use of the \code{\link{gnlr}} function, the vector of stimuli is defined as a global variable. It's recommended to avoid global variables for better code maintainability and to prevent potential conflicts. This implementation will be changed in future releases.
#' When guess and lapse parameters are present, the fitted coefficients returned by the \code{glnr} function undergo a log-transformation. This transformation is applied as a result of fitting the parameters of the model using the exponential function, which ensures that the coefficients remain positive.
#' 
#'
PsychFunction_gnlm <- function(model_glm, response, stimuli, link, data, guess, lapse){
  
  stopifnot(is.character(response), is.character(stimuli))
  
  response <- parse(text = paste("cbind(", response[1], ",", response[2], ")"))
  
  glm_coeff <- summary(model_glm)$coefficients[,1]
  start_estimate <- c(-glm_coeff[1]/glm_coeff[2], ifelse(link == "weibull", 2, 1/glm_coeff[2]))
  
  lambda <- if (isTRUE(lapse)) runif(1, min = 0, max = 0.05) else if (is.numeric(lapse)) lapse else FALSE
  gamma <- if (isTRUE(guess)) runif(1, min = 0, max = 0.05) else if (is.numeric(guess)) guess else FALSE
  
  if (is.numeric(gamma) && !is.numeric(lambda)){
    start_estimate <- c(start_estimate, log(gamma))
  }else if (!is.numeric(gamma) && is.numeric(lambda)){
    start_estimate <- c(start_estimate, log(lambda))
  }else if(is.numeric(gamma) && is.numeric(lambda)){
    start_estimate <- c(start_estimate, log(gamma), log(lambda))
  }

  #make sure you don't overwrite existing variables
  if (exists("x_values", envir = .GlobalEnv)) {
    if (identical(parent.frame(2), .GlobalEnv)){
      decide <- readline("The variable x_values in the global environment will be overwritten and deleted. Do you wish to continue? (yes/no): ")
    }else{
      decide <- "y"
    }
    if (tolower(decide) %in% c("yes", "y")) {
      rmGlobalVar()
    }else{
      stop("Function terminated by user.")
    }
  } else {
    message("x_values will be created in your global environment and quickly deleted - just so you know.")
  }
  source(here("R", "global.R"))
  setGlobalVar(data, stimuli) #x_values defined as global variable due to gnlm syntax. 
  switch_function <- switch_mu_function(func_name = link, gamma, lambda)
  model_gnlm <- gnlr(y = with(data, eval(response)), distribution = "binomial",
                     mu = switch_function, pmu = start_estimate)
  
  rmGlobalVar()
  
  gnlr_coeff <- model_gnlm$coefficients
  
  if (length(gnlr_coeff) > 2) {
    gnlr_coeff[3:length(gnlr_coeff)] <- log(gnlr_coeff[3:length(gnlr_coeff)])
  }
  
  return(list(model = model_gnlm, gnlr_coeff = gnlr_coeff))
}

#' Internal Function: Switch mu Function
#'
#' This function switches between different mu functions based on the provided parameters.
#'
#' @param func_name A string specifying the function being fitted. Possible options are 'probit' for the cumulative normal distribution, 'logit' for the cumulative logit distribution, 'weibull' for cumulative Weibull distribution. 
#' @param gamma,lambda parameters indicating whether to include guessing and lapse parameters, respectively. If parameters are FALSE, they are not included in the model. If numeric, they are the starting estimates used in \code{\link{gnlr}}.
#' 
#' @note
#' Instead of directly fitting guess and lapse parameters, this function uses \code{exp(p[3])} and \code{exp(p[4])}. By doing so, we force the estimations to be non-negative. 
#' 
#' 
#' @seealso 
#' [\code{\link{PsychFunction}}] [\code{\link{PsychFunction_gnlm}}] [\code{\link{gnlr}}]
#'
#' @return A function representing the selected mu function.
#' @importFrom stats pweibull
#' 
switch_mu_function <- function(func_name, gamma, lambda) {
  #mu <- function(p) pnorm(x_values, mean = p[1], sd = p[2])
  if (isFALSE(gamma) && isFALSE(lambda)){
    switch(func_name,
          probit = function(p) pnorm(x_values, mean = p[1], sd = p[2]),
          logit = function(p) plogis(x_values, location = p[1], scale = p[2]),
          weibull = function(p) pweibull(x_values, scale = p[1], shape = p[2]))
  }else if (is.numeric(gamma) && isFALSE(lambda)){
    switch(func_name,
           probit = function(p) exp(p[3]) + (1 - exp(p[3])) * pnorm(x_values, mean = p[1], sd = p[2]),
           logit = function(p) exp(p[3]) + (1 - exp(p[3])) * plogis(x_values, location = p[1], scale = p[2]),
           weibull = function(p) exp(p[3]) + (1 - exp(p[3])) * pweibull(x_values, scale = p[1], shape = p[2]))
  }else if (isFALSE(gamma) && is.numeric(lambda)){
    switch(func_name,
           probit = function(p) (1 - exp(p[3])) * pnorm(x_values, mean = p[1], sd = p[2]),
           logit = function(p) (1 - exp(p[3])) * plogis(x_values, location = p[1], scale = p[2]),
           weibull = function(p) (1 - exp(p[3])) * pweibull(x_values, scale = p[1], shape = p[2]))
  }else if(is.numeric(gamma) && is.numeric(lambda)){
    switch(func_name,
           probit = function(p) exp(p[3]) + (1 - exp(p[3]) - exp(p[4])) * pnorm(x_values, mean = p[1], sd = p[2]),
           logit = function(p) exp(p[3]) + (1 - exp(p[3]) - exp(p[4])) * plogis(x_values, location = p[1], scale = p[2]),
           weibull = function(p) exp(p[3]) + (1 - exp(p[3]) - exp(p[4])) * pweibull(x_values, scale = p[1], shape = p[2]))
  }
}

#' Internal Function: Set global variable
#' 
#' Setting a global variable is necessary in order to correctly define the function to fit with gnlm. 
#' 
#' @param data dataset
#' @param stimuli A string of the stumuli variable.
#' 
#' @seealso [PsychFunction_gnlm()]
#' @noRd
#' @return Invisible NULL.
#'
#' @keywords internal

setGlobalVar <- function(data, stimuli) {
  unlockBinding("x_values", globalenv())
  #globalenv()$x_values <- data[, stimuli]
  assign("x_values", data[, stimuli], envir = .GlobalEnv)
  #x_values <<- data[, stimuli]
  # Lock the binding after modification
  lockBinding("x_values", globalenv())
  #invisible(NULL)
}

#' Internal Function: Delete global variable
#' 
#' Remove the global variable
#' @noRd
#' 
rmGlobalVar <- function(){
  rm("x_values", envir = .GlobalEnv)
}