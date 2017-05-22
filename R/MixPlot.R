#' Plotting univariable GLMM
#'
#' A function to plot binomial data and the fitted GLMM (object of class xplode).
#' 
#'  @param xplode.obj an object of class xplode
#'  @param pf for multivariable GLMM including one factorial predictor,
#'   the level number to be plotted 
#'  @param p05line logical, should be an horizontal and a vertical line added? 
#'   the horizontal line is fixed at P(Y = 1) = 0.5.
#'  @param x.ref if p05line = T, this is the position of the vertical line on the x axis
#'  @param x.range a vector of length two specifying the range for model predictions
#'  @param col logical, if TRUE a different color will be used for different clusters/participants
#'  @param x.label,y.label label for the x and the y axes
#'   
#'  @details The function is currently only working with GLMM including maximum three 
#'  random effects (random intercept, random slope and covariance of the two)
#'  
#'  @return a data.frame object including the intercept and slope of each participant
#'   (algebraic sum of the fixed effects and the modes of the random effects) and the 
#'   color number for the plot. 
#'
#'  @examples
#' datafr.1 <- PsySimulate(fixeff = c(-7.5, 0.0875), nsubject = 6, constant = T)                               
#' formula.mod = cbind(Longer, Total - Longer) ~ X + (1 + X| Subject)
#' mod1 <- glmer(formula = formula.mod, family = binomial(link = "probit"), data = datafr.1)
#' define.mod = list(pf1 = list(intercept = 1, slope = 2))
#' xplode.mod1 = xplode.mer(model = mod1, name.cont = "X", define.pf = define.mod)
#' myplot = MixPlot(xplode.mod1, pf = 1,  p05line = F, x.ref = 60, x.range = c(0,160),
#'  col = T,x.label = "Stimulus Intensity", y.label = "Predicted Response")
#'
#' @export
MixPlot = function(xplode.obj, pf = 1,  p05line = F, x.ref,
                   x.range, col = F,
                   x.label = "Stimulus Intensity", y.label = "Predicted Response"){
  
  #numebr of subjects
  nsubjects = xplode.obj$Groups
  
  #only one condition
  if(xplode.obj$n.psych.fun > 1){
    newdata = subset.data.frame(xplode.obj$model.frame,
                                xplode.obj$model.frame[ ,xplode.obj$factor.col] == xplode.obj$factor.levels[pf])
  }else{
    newdata = xplode.obj$model.frame
  }
  newdata$size = rowSums(newdata[,1])
  
  #colors and pch 
  if(col == T){
    if(length(palette()) < nsubjects){palette(rainbow(nsubjects))}		
    newdata$pointcolors = newdata[,xplode.obj$Groups.colnames]
    curvecolors = 1:nsubjects
  }else{
    newdata$pointcolors = 1
    curvecolors = rep(1, nsubjects)
  }

  #BLUPS (random modes)
  fixef.pf = cbind(xplode.obj$psychometrics[[pf]]$intercept["Estimate"],
                   xplode.obj$psychometrics[[pf]]$slope["Estimate"])

  if(xplode.obj$multirand == FALSE){
    estimates = data.frame(matrix(numeric(), nrow = nrow(ranef(model)[[1]]), ncol = 3))
    estimates[,1] = xplode.obj$ranef + fixef.pf[1]
    estimates[,2] = rep(fixef.pf[2], times = nrow(estimates))
    estimates[,3] = curvecolors
  }else{
    estimates = t(apply(X = xplode.obj$ranef[[1]], MARGIN = 1, FUN = function(X){X + fixef.pf}))
    estimates = data.frame(cbind(estimates, curvecolors))
  }
  
  #see xplode.mer$moddel.frame
  plot(y =  newdata[,1][,1]/newdata$size,
       x = newdata[,xplode.obj$cont.col],  xlim = x.range,
       xlab = x.label, ylim = c(0,1), ylab = y.label,
       col = newdata$pointcolors)		 													
  
  if (p05line == T){
    abline(h = 0.5, col = "lightgrey", lty = "dashed")
    abline(v = x.ref, col = "lightgrey", lty = "dashed")
  }
  
  curve.probit(estimates, x.range[1], x.range[2])	
  
  palette("default")	
  return(estimates)
}