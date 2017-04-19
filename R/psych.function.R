#psych.function fit a glm to binary data and plot the fitted function
#AM did it on 25.07.2013
#AM modified on 30.07.2015: PSE and JND in the outcome
#To do: multivariable model plotting

#ps.formula		formula object such as as.formula("cbind(away, towards) ~ velocity")
#ps.link		character. link function
#ps.data  	data.frame including the binary response variables and predictor
#x.range		the range over which the function will be plotted.
#ps.col     character. Color of the line
#ps.lines   logical, indicatinf if a line of the model fit should be added to the existing plot 

psych.function = function(ps.formula, ps.link, ps.data, x.range = c(NA, NA), ps.x = NA, 
                          ps.lines = F, ps.col = "black", br = F){
	myfit = vector("list", 2)
  ps.terms = terms(ps.formula)
  
  model.glm = glm(formula = ps.formula, family = binomial(link = ps.link), data = ps.data)
  brflag = ifelse(trunc(max(model.glm$fitted.values)) == 1 | trunc(min(model.glm$fitted.values)) == 0, T, F)
  myfit[[3]] = brflag
  
  if(br == T & brflag == T){
    require(brglm)
    model.glm = brglm(formula = ps.formula, family = binomial(link = ps.link), data = ps.data)
    warning("Binomial-response GLMs fitted using the bias-reduction method (brglm)")
  }
	
  myfit[[1]] = model.glm
  
  if(ps.link == "probit"){
    myfit[[2]] = delta.psy.probit(model.glm)
    if(ps.lines == T){
      segments(x0 = myfit[[2]][1,3], y0 = 0.5, x1 = myfit[[2]][1,4], y1 = 0.5, col = ps.col)
    }
  }else{
    myfit[[2]] = NA
    warning("Use the probit link function to get the estimate of the PSE and the JND")
  }
  
	if(ps.lines == T){
    if(is.na(ps.x)){
      ps.x <- data.frame(seq(x.range[1], x.range[2], length.out = 1000))
      names(ps.x) = attr(ps.terms, "term.labels")
    }
		lines(x = seq(x.range[1], x.range[2], length.out = 1000),
		 y = predict(object = model.glm, newdata = data.frame(ps.x), type = "response"),
		 col = ps.col)
		}
	return(myfit)
	}
