delta.psy.probit <- function(model, alpha = 0.05, lme4 = F){

  if(lme4 == F){
  pse <- -model$coef[1]/model$coef[2]
  BETA <- model$coef[2]} else {
  #		if extracting form mer, check lme4 version!	
  fixed.par = getME(model, "beta")	
  pse <- -(fixed.par[1]/fixed.par[2])
  BETA <- fixed.par[2]}
  
  #		if extracting form mer, check lme4 version!
  var.alpha <- vcov(model)[1,1]
  var.beta <- vcov(model)[2,2]
  cov.alpha.beta <- vcov(model)[2,1]
  
  var.pse <- (1/BETA^2)*(var.alpha + (2*pse*cov.alpha.beta)+(pse^2*var.beta))   #PSE
  inferior.pse <- pse - (qnorm(1 - (alpha/2))*sqrt(var.pse))
  superior.pse <- pse + (qnorm(1 - (alpha/2))*sqrt(var.pse))

  jnd <- qnorm(0.75) * (1/BETA)
  var.jnd <- (qnorm(0.75) * (-1/BETA^2))^2 * var.beta                           #JND
  inferior.jnd <- jnd - (qnorm(1 - (alpha/2))*sqrt(var.jnd))
  superior.jnd <- jnd + (qnorm(1 - (alpha/2))*sqrt(var.jnd))

  output <- matrix(rbind(c(pse, sqrt(var.pse), inferior.pse, superior.pse),
    c(jnd, sqrt(var.jnd), inferior.jnd, superior.jnd)), nrow = 2,
    dimnames = list(param <- c("pse", "jnd"), statistics <- c("Estimate",
    "Std. Error", "Inferior", "Superior")))

return(output)}
