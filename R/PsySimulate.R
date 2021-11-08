#' Simulate psychophysical data
#'
#' The function simulates data of a typical psychophysics experiment. For each simulated 
#' participant, the function returns the following information: individual slope and intercept 
#' coefficients, given the fixed and random effects parameters provided as input; summary of 
#' the simulated binomial response to a range of intensity levels between a specified range.
#'
#' @param fixeff  array of fixed effects. First item is the intercept, second
#' item is the slope.
#' @param raneff array of random effects. First item is the intercept, second
#' item is the covariance, third item is the slope.
#' @param nsubjects number of subjects to simulate. Default is 8.
#' @param ntrials number of trials for each stimulus level. Default is 40.
#' @param nintervals number of stimulus levels. Default is 9.
#' @param xint range of the stimulus interval. Default is c(40,120)
#' @param constant logical. If set to FALSE, stimulus levels are randomly generated,
#' uniformly distributed values within the selected interval.
#' otherwise, the X interval is divided in  intervals of constant
#' length. Default is TRUE.
#'
#' @return \code{PsySimulate} returns a simulated dataset. If no input arguments are specified, the function returns 
#' a dataset with the same characteristics as \code{\link{simul_data}}.
#'
#' @seealso \code{\link{PsychShape}} for plotting a psychometric function given PSE and JND. 
#'
#' @examples
#' datafr.S1 <- PsySimulate(fixeff = c(0, 1), xint = c(-5,5), 
#' nsubject = 1, ntrials = 60, nintervals = 10, constant = FALSE)
#' library(ggplot2)
#' g <- ggplot(datafr.S1, aes(X,Longer/Total)) + geom_point()
#' 
#' PsychShape(pse = 0, jnd = qnorm(0.75)/1, ps.link = "probit", 
#' x.range = c(-5,5), addTo = g, ps.color = "red")
#'
#' @importFrom mnormt rmnorm
#' @importFrom Matrix nearPD
#' @export
#'
PsySimulate <- function(fixeff = c(-7, 0.0875), raneff = c(2.4, -0.002, 2e-06), nsubjects = 8, ntrials = 40, nintervals = 9, xint = c(40, 120), constant = T) {
    
    
    if (constant == T) {
        xval = seq(from = min(xint), to = max(xint), by = (max(xint) - min(xint))/(nintervals - 1))
    } else {
        xval = sort(runif(n = nintervals, min = min(xint), max = max(xint)))
    }
    
    
    if (nsubjects == 1) {
        GRID <- data.frame(X = xval, Intercept = fixeff[1], Slope = fixeff[2], Total = ntrials, Subject = 1)
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
        
        GRID <- data.frame(lapply(buffer,unlist), Total = ntrials)
        GRID <- GRID[order(GRID$Subject, GRID$X),]
    }
    
    GRID$prob <- apply(GRID, MARGIN = 1, FUN = function(X) {
        pnorm(q = X["X"], mean = -X["Intercept"]/X["Slope"], sd = abs(1/X["Slope"]))
    })
    
    
    GRID$Longer <- apply(GRID, MARGIN = 1, FUN = function(X) rbinom(n = 1, prob = ifelse(X["Slope"] > 0, X["prob"], 1-X["prob"]), size = X["Total"]))
    
    
    keep <- c("X", "Intercept", "Slope", "Longer", "Total", "Subject")
    datafr <- GRID[keep]
    
    return(datafr)
    
}
