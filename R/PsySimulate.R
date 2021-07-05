#' Simulate psychophysical data
#'
#' Given the arrays of fixed and random effects, as well as the covariance,
#' and the characteristic of the simulated experiment (i.e., ) the function
#' simulates a dataset in which for each subject the following information is
#' provided: the slope and intercept value of the psychometric function, and
#' the simulated responses to the stimulus levels that fit that function.
#'
#' @param fixeff  Array of fixed effects. First item is the intercept, second
#' item is the slope.
#' @param raneff Array of random effects. First item is the intercept, second
#' item is the covariance, third item is the slope.
#' @param nsubjects Number of subjects to simulate data for. Default is 8.
#' @param ntrials Number of trials for each stimulus level. Default is 40.
#' @param nintervals Number of stimulus levels. Default is 9.
#' @param xint Range of the stimulus interval. Default is c(40,120)
#' @param constant If set to FALSE, stimulus levels are randomly generated,
#' uniformly distributed values within the selected interval.
#' If constant = TRUE, the X interval is divided in  intervals of constant
#' length. Default is TRUE.
#'
#' @return The simulated dataset
#'
#' @examples
#' #simulate dataset (one subject)
#' datafr.S1 <- PsySimulate(nsubject = 1, constant = TRUE)
#'
#' @importFrom mnormt rmnorm
#' @importFrom Matrix nearPD
#' @importFrom stats rbinom runif
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
