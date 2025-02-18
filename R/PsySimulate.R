#' Simulate psychophysical data
#'
#' The function simulates data from a typical psychophysics experiment using a 2-alternative forced choice task and the method of constant stimuli. 
#' For each simulated participant, the function returns the following information: individual slope and intercept coefficients, based on the fixed and random effects parameters provided as input; summary of the simulated binomial response across a range of intensity levels within a specified stimulus range; individual lapse and guess rates, if applicable.
#'
#' @param fixeff  Numeric array of fixed effects. The first item is the intercept, the second element is the slope.
#' @param raneff Numeric array of random effects. The first element is the intercept, the second is the covariance, the third is the slope variance.
#' @param nsubjects Integer. Number of subjects to simulate. Default is 8.
#' @param ntrials Integer. Number of trials for each stimulus level. Default is 40.
#' @param nintervals Integer. Number of stimulus levels. Default is 9.
#' @param xint Numeric array specifying the range of stimulus intensity. Default is c(40,120)
#' @param constant Logical. If TRUE (defualt), stimulus levels are evenly spaced across `xint`. If FALSE, stimulus levels are randomly generated within the interval.
#' @param lapse Logical or numeric. If FALSE (default), no lapse rate is applied. If TRUE, a random lapse rate is drawn from `lapse_range`. If a numeric value is provided, all subjects will have the same lapse rate.
#' @param guess Logical or numeric. If FALSE (default), no guessing rate is applied. If TRUE, a random guessing rate is drawn from `guess_range`. If a numeric value is provided, all subjects will have the same guess rate.
#' @param lapse_range Numeric array defining the minimum and maximum lapse rates when `lapse = TRUE`. Default is c(0, 0.05).
#' @param guess_range Numeric array defining the minimum and maximum guessing rates when `guess = TRUE`. Default is c(0.05, 0.10).
#' 
#' @return A data frame containing simulated psychophysical data with the following columns:
#' \itemize{
#'   \item \code{Subject} - Subject identifier.
#'   \item \code{X} - Stimulus intensity levels.
#'   \item \code{Intercept} - Individual intercept values.
#'   \item \code{Slope} - Individual slope values.
#'   \item \code{Gamma} (optional) - Guess rate, included if `guess` is not FALSE.
#'   \item \code{Lambda} (optional) - Lapse rate, included if `lapse` is not FALSE.
#'   \item \code{Longer} - Number of "Longer" responses at each stimulus level.
#'   \item \code{Total} - Total number of trials per stimulus level.
#' }
#'
#' @seealso \code{\link{PsychShape}} for plotting a psychometric function given PSE and JND. \code{\link{simul_data}} for a dataset simulated with the function.
#'
#' @examples
#' datafr.S1 <- PsySimulate(fixeff = c(0, 1), xint = c(-5,5), 
#' nsubject = 1, ntrials = 60, nintervals = 10, constant = FALSE)
#' 
#' simul_data <- PsySimulate(ntrials = 160, nsubjects = 10, guess = TRUE, lapse = TRUE)
#'
#' @importFrom mnormt rmnorm
#' @importFrom Matrix nearPD
#' @importFrom stats runif pnorm rbinom 
#' @export
#'
PsySimulate <- function(fixeff = c(-7, 0.0875), raneff = c(2.4, -0.002, 2e-06), nsubjects = 8, ntrials = 40, nintervals = 9, xint = c(40, 120), constant = T, lapse = FALSE, guess = FALSE, lapse_range = c(0,0.05), guess_range = c(0, 0.05)) {
  
  
  if (nsubjects > 1 && is.numeric(lapse) && lapse > 0) {
    warning("Setting lapse rate as a number > 0 with more than one subject. Each subject will have the same lapse rate.")
  }
  
  if (nsubjects > 1 && is.numeric(guess) && guess > 0) {
    warning("Setting guess rate as a number > 0 with more than one subject. Each subject will have the same guess rate.")
  }
  
  
  if (constant) {
    xval = seq(from = min(xint), to = max(xint), by = (max(xint) - min(xint))/(nintervals - 1))
  } else {
    xval = sort(runif(n = nintervals, min = min(xint), max = max(xint)))
  }
  
  
  if (nsubjects == 1) {
    Lambda <- if (!lapse) 0 else if (isTRUE(lapse)) runif(1, min = lapse_range[1], max = lapse_range[2]) else lapse
    Gamma <- if (!guess) 0 else if (isTRUE(guess)) runif(1, min = guess_range[1], max = guess_range[2]) else guess
    GRID <- data.frame(X = xval, Intercept = fixeff[1], Slope = fixeff[2], Total = ntrials, Subject = 1, Lambda = Lambda, Gamma = Gamma)
  }else{
    #Variance and covariance matrix of the random effects
    vcov.matrix = matrix(raneff[c(1, 2, 2, 3)], nrow = 2)
    # nearPD: near positive defined matrix
    mer.vcov = nearPD(vcov.matrix)$mat  
    
    # independent samples from multivariate normal distribution
    buffer <- rmnorm(n = nsubjects, mean = rep(0, 2), varcov = mer.vcov)
    colnames(buffer) <- c("sample.interc", "sample.slope")
    
    buffer <- data.frame(buffer, Subject = seq(1:nsubjects))
    
    buffer$Intercept <- buffer$sample.interc + fixeff[1]
    buffer$Slope <- buffer$sample.slope + fixeff[2]
    buffer$X <- list(xval)
    buffer$Lambda <- if (!lapse) 0 else if (isTRUE(lapse)) runif(nsubjects, min = lapse_range[1], max = lapse_range[2]) else lapse
    buffer$Gamma <- if (!guess) 0 else if (isTRUE(guess)) runif(nsubjects, min = guess_range[1], max = guess_range[2]) else guess
    
    
    GRID <- data.frame(lapply(buffer,unlist), Total = ntrials)
    GRID <- GRID[order(GRID$Subject, GRID$X),]
  }
  
  GRID$prob <- apply(GRID, MARGIN = 1, FUN = function(X) {
    X["Gamma"] + (1 - X["Gamma"] - X["Lambda"]) * pnorm(q = X["X"], mean = -X["Intercept"]/X["Slope"], sd = abs(1/X["Slope"]))
  })
  
  GRID$Longer <- apply(GRID, MARGIN = 1, FUN = function(X) rbinom(n = 1, prob = ifelse(X["Slope"] > 0, X["prob"], 1-X["prob"]), size = X["Total"]))
  
  
  keep <- c("Subject", "X", "Intercept", "Slope", "Longer", "Total")
  
  if (guess != FALSE) {
    keep <- append(keep, "Gamma", after = match("Slope", keep))
  }
  if (lapse != FALSE) {
    keep <- append(keep, "Lambda", after = match("Gamma", keep, nomatch = match("Slope", keep)))
  }
  
  datafr <- GRID[keep]
  
  return(datafr)
  
}
