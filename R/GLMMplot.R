#' @importFrom grDevices palette rainbow
#' @importFrom graphics abline plot
#' 
GLMMplot = function(dataframe, X.col, Yes.col, Total.col, Subject.col, estimates, lme4 = F, model,
   		x.label = "Stimulus Intensity", y.label = "Predicted Response",
   		p05line = F, x.ref = mean(dataframe[,X.col]),
   	 	x.from = min(dataframe[,X.col]), x.to = max(dataframe[,X.col]), col = F){
	
	#numebr of subjects
	nsubjects = nlevels(factor(dataframe[, Subject.col]))
	
	#colors and pch		
	if(col == T){
		if(length(palette()) < nsubjects){palette(rainbow(nsubjects))}		
		curvecolors = 1:nsubjects
		dataframe$pointcolors = dataframe[, Subject.col]}else{
				 	curvecolors = rep(1, nsubjects)
				 	dataframe$pointcolors = rep(1, times = nrow(dataframe))}
	dataframe$pch = rep(1, times = nrow(dataframe))
	
	#for plotting beta > 0 and beta < 0
	#----------------------------------------------------------------------------#	
	plot.cdf = function(X){
		BETAplus = which(X[,2] > 0) 
		Xplus = X[BETAplus,]
		Xminus = X[-BETAplus,]
		x = NA
		
		if(nrow(Xplus) > 0){	
		apply(X = Xplus, MARGIN = 1,
		 FUN = function(X) {curve(expr = pnorm(x, mean = -X[1]/X[2], sd = 1/X[2] ),
		 	 from = x.from, to = x.to,
		 	  #col = "black", 
		 	  col = X[3],
		 	  add = T)}
		 )}
		 
		if(nrow(Xminus) > 0){
		apply(X = Xminus, MARGIN = 1,
		 FUN = function(X) {curve(expr = 1 - pnorm(x, mean = -X[1]/X[2], sd = 1/abs(X[2] )),
		 	from = x.from, to = x.to,
		 		#col = "red", 
		 	col = X[3],
		 	add = T)}
		 )}
		return(BETAplus)
		}
	#----------------------------------------------------------------------------#
	
	#BLUPS (random modes)
	if(lme4 == T){
		if(dim(ranef(model)[[1]])[2] == 1){
			estimates = data.frame(matrix(numeric(), nrow = nrow(ranef(model)[[1]]), ncol = 3))
			estimates[,1] = ranef(model)[[1]] + fixef(model)[1]
			estimates[,2] = rep(fixef(model)[2], times = nrow(estimates))
			estimates[,3] = curvecolors}else{
			estimates = t(apply(X = ranef(model)[[1]], MARGIN = 1, FUN = function(X){X + fixef(model)}))
			estimates = cbind(estimates, curvecolors)}
			}

	
	plot(y = dataframe[,Yes.col]/dataframe[,Total.col], x = dataframe[,X.col],  xlim =c(x.from, x.to),
	 xlab = x.label, ylim = c(0,1), ylab = y.label,
	 col = dataframe$pointcolors, pch = dataframe$pch)		 													
	
	if (p05line == T){
		abline(h = 0.5, col = "lightgrey", lty = "dashed")
		abline(v = x.ref, col = "lightgrey", lty = "dashed")}
		
	plot.cdf(estimates)	
	
	palette("default")	
	return(list(dataframe, estimates))
}

#' Plotting univariable GLMM
#'
#' Plot binomial data and the fitted GLMM (object of class xplode).
#' 
#' @param xplode.obj an object of class xplode
#' @param pf integer: for multivariable GLMM including one factorial predictor,
#' the level number to be plotted 
#' @param p05line logical, should be an horizontal and a vertical line added? 
#' the horizontal line is fixed at P(Y = 1) = 0.5.
#' @param x.ref if p05line = T, this is the position of the vertical line on the x axis
#' @param x.range a vector of length two specifying the range for model predictions
#' @param col logical, if TRUE a different color will be used for different clusters/participants
#' @param x.label,y.label label for the x and the y axes. If not specified, x.labels = ""Stimulus Intensity", 
#' y.label = "Predicted Response"
#'   
#' @note The function is currently only working with GLMM including maximum three 
#' random effects (random intercept, random slope and covariance of the two)
#'  
#' @return a data.frame object including the intercept and slope for each participant
#' (algebraic sum of the fixed effects and the modes of the random effects) and the 
#' color number for the plot. 
#' 
#' @seealso \code{\link{xplode}} objects of class \code{xplode.obj}.
#' 
#' @keywords DeltaMethod GLMM Plotting
#'
#' @examples
#' library(lme4)
#' data(vibro_exp3)
#' formula.mod <- cbind(faster, slower) ~ speed + (1 + speed| subject)
#' mod <- glmer(formula = formula.mod, family = binomial(link = "probit"),
#'               data = vibro_exp3[vibro_exp3$vibration == 0,])
#' define.mod <- list(pf1 = list(intercept = 1, slope = 2))
#' xplode.mod <- xplode(model = mod, name.cont = "speed", define.pf = define.mod)
#' myplot <- MixPlot_1.0(xplode.mod, pf = 1,  p05line = FALSE, x.ref = 8.5, x.range = c(1,16),
#'                   col = TRUE, x.label = "Stimulus Speed", y.label = "Predicted Response")
#'
#' @export
#' 
MixPlot_1.0 <- function(xplode.obj, pf = 1, p05line = F, x.range, x.ref,
                    col = F, x.label = "Stimulus Intensity", y.label = "Predicted Response"){
  
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
    estimates = data.frame(matrix(numeric(), nrow = nrow(xplode.obj$ranef[[1]]), ncol = 3))
    estimates[,1] = xplode.obj$ranef[[1]] + fixef.pf[1]
    estimates[,2] = rep(fixef.pf[2], times = nrow(estimates))
    estimates[,3] = curvecolors
  }else{
    estimates = t(apply(X = xplode.obj$ranef[[1]], MARGIN = 1, FUN = function(X){X + fixef.pf}))
    estimates = data.frame(cbind(estimates, curvecolors))
  }
  
  #see xplode$moddel.frame
  plot(y =  newdata[,1][,1]/newdata$size,
       x = newdata[,xplode.obj$cont.col],  xlim = x.range,
       xlab = x.label, ylim = c(0,1), ylab = y.label,
       col = newdata$pointcolors)		 													
  
  if (p05line == T){
    abline(h = 0.5, col = "lightgrey", lty = "dashed")
    abline(v = x.ref, col = "lightgrey", lty = "dashed")
  }
  
  CurveProbit(estimates, x.range[1], x.range[2])	
  
  palette("default")	
  return(estimates)
}
