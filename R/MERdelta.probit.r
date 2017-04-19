MERdelta.probit <- function(xplode.obj, alpha = 0.05){
  
  #check if link = probit
  if(xplode.obj$family$link != "probit"){
      output = NA
      print("Use a probit link function")
  }else{
      n.pf = length(xplode.obj$psychometrics)
      output = vector("list", length = n.pf)
      names(output) = names(xplode.obj$psychometrics)
      
      for(i in 1:n.pf){
        #copy all the variables in temporary objects
        pse <- -(xplode.obj$psychometrics[[i]]$intercept[1]/
                                      xplode.obj$psychometrics[[i]]$slope[1])
        slope <- xplode.obj$psychometrics[[i]]$slope[1]
        
        var.intercept <- xplode.obj$psychometrics[[i]]$intercept[2]
        var.slope <- xplode.obj$psychometrics[[i]]$slope[2]
        
        #cov(alpha, slope): for all pfs, is approximated to the cov(alpha1, slope1)
        cov.intercept.slope <- xplode.obj$psychometrics$pf1$cov
        
        #compute all the other variables
        var.pse <- (1/slope^2)*(var.intercept + (2*pse*cov.intercept.slope)+(pse^2*var.slope))   #PSE
        inferior.pse <- pse - (qnorm(1 - (alpha/2))*sqrt(var.pse))
        superior.pse <- pse + (qnorm(1 - (alpha/2))*sqrt(var.pse))
        
        jnd <- qnorm(0.75) * (1/slope)
        var.jnd <- (qnorm(0.75) * (-1/slope^2))^2 * var.slope                           #JND
        inferior.jnd <- jnd - (qnorm(1 - (alpha/2))*sqrt(var.jnd))
        superior.jnd <- jnd + (qnorm(1 - (alpha/2))*sqrt(var.jnd))
        
        output[[i]] <- matrix(rbind(
                                c(pse, sqrt(var.pse), inferior.pse, superior.pse),
                                c(jnd, sqrt(var.jnd), inferior.jnd, superior.jnd)),
                              nrow = 2,
                              dimnames = list(param <- c("pse", "jnd"),
                                              statistics = c("Estimate","Std. Error", "Inferior", "Superior")))
      }
    }
    

  return(output)
}
