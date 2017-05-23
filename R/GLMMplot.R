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