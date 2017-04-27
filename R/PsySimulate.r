PsySimulate = function(fixeff = c(-7, 0.0875), raneff = c(2.4, -20e-04, 2e-06),
                         nsubjects = 8, pps = 9, ntrials = 40,
                         xint = c(40, 120), constant = F){
  
  
  if(constant == T){xval = seq(from = min(xint), to = max(xint), by = (max(xint)-min(xint))/(pps-1))
  }else{xval = sort(runif(n = pps, min = min(xint), max = max(xint)))}
  
  
  if(nsubjects == 1){
    
    r.interc <- rep(fixeff[1], pps)
    r.slope <- rep(fixeff[2], pps)
    
    
    GRID <- data.frame(cbind(xval, r.interc, r.slope))
    names(GRID) <- c("X", "alpha", "beta")
    
    probab <- apply(GRID, MARGIN = 1,
                    FUN = function(X){pnorm(q = X[1], mean = -X[2]/X[3], sd = 1/X[3])})
    probab <- data.frame(rep(ntrials, pps), probab)
    names(probab) <- c("size", "prob")
    
    
    response <- apply(probab, MARGIN = 1, FUN = function(X) rbinom(n = 1,
                                                                   prob = X[2], size = X[1]))
    
    datafr <- data.frame(cbind(GRID[,1:3], response, probab[,1],
                               as.factor(rep(1, pps))))     
    
    names(datafr) <- c("X", "Intercept", "Slope", "Longer", "Total","Subject")
    return(datafr)
    
  }else{
    #library(Matrix)
    #library(mnormt)
    vcov.matrix = matrix(raneff[c(1,2,2,3)], nrow = 2)	
    mer.vcov = nearPD(vcov.matrix)	 							#Variance and covariance matrix of the random effects 
    #nearPD: near positive defined matrix
    buffer <- matrix(numeric(), nrow = nsubjects, ncol = 2)
    
    for(i in 1:nsubjects){buffer[i,1:2] <- rmnorm(n = 1, 							#independent samples from multivariate normal distribution
                                                  mean = rep(0, 2), varcov = as.matrix(mer.vcov$mat))} 
    
    
    row.dataframe = pps*nsubjects									#The number of rows of the dataframe
    r.interc <- matrix(data = numeric(), nrow = row.dataframe,
                       ncol = 5, dimnames = list(c(NULL), c("Subj", "Mean", "Str.Dev", "Sample.Val",
                                                            "Fix.Eff")) )
    r.interc[,1] <- sort(rep(1: nsubjects, pps))
    r.interc[,2] <- rep(0, row.dataframe)
    r.interc[,3] <- rep(raneff[1], row.dataframe)
    r.interc[,4]<- rep(buffer[,1], each = pps)
    r.interc[,5] <- rep(fixeff[1], row.dataframe)
    
    r.slope <- matrix(data = numeric(), nrow = row.dataframe,
                      ncol = 5, dimnames = list(c(NULL), c("Subj", "Mean", "Str.Dev", "Sample.Val",
                                                           "Fix.Eff")) )
    r.slope[,1] <- sort(rep(1: nsubjects, pps))
    r.slope[,2] <- rep(0, row.dataframe)
    r.slope[,3] <- rep(raneff[3], row.dataframe)									#random slope -> raneff[2] is covariance
    r.slope[,4] <- rep(buffer[,2], each = pps)
    r.slope[,5] <- rep(fixeff[2], row.dataframe)
    
    
    xval2 = rep(xval, nsubjects)
    GRID <- data.frame(cbind(xval2,
                             r.interc[,4] + r.interc[,5], r.slope[,4] + r.slope[,5]))
    names(GRID) <- c("X", "alpha", "beta")                           
    
    
    probab <- apply(GRID, MARGIN = 1,
                    FUN = function(X){pnorm(q = X[1], mean = -X[2]/X[3], sd = 1/X[3])})
    probab <- data.frame(rep(ntrials, row.dataframe), probab)
    names(probab) <- c("size", "prob")
    
    
    response <- apply(probab, MARGIN = 1, FUN = function(X) rbinom(n = 1,
                                                                   prob = X[2], size = X[1]))
    
    
    datafr <- data.frame(cbind(GRID[,1:3], response, probab[,1],
                               as.factor(r.interc[,1])))                             
    names(datafr) <- c("X", "Intercept", "Slope", "Longer", "Total","Subject")
    
    return(datafr)
  }
}