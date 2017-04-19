#xplode.mer extract values and parametrs from a mer object. The FUN defines the statistics of interests
#To DO: adjust at line 115 for binary (single col) data
xplode.mer = function(model, name.cont = NA, name.factor = NA, names.response = NA,
                      define.pf = list(pf1 = list(intercept = 1, slope = 2)))
  {
	xplode = vector("list", 25)
	names(xplode) = c( "fixef", "fixef.vcov", "n.psych.fun",
                     "factor.col", "factor.colname", "factor.parnames",
                      "cont.col", "cont.colname",
                      "psychometrics", "ranef.stddev", "ranef.VarCov", "multirand", 
	                    "Groups","Gp","Groups.levels", "Groups.colnames", "flist", "nobs",
                      "size", "response.colnames", "model.frame", "model.matrix", "rps",
                      "formula", "family")
  
  
  #1)fixed effects--------------------------
	xplode$fixef = fixef(model)
	
  #Variance-Covariance Matrix of the fixed effects
	xplode$fixef.vcov = vcov(model)
	
  #Number of psychometric functions
	xplode$n.psych.fun = length(define.pf)
  
	#Number of psychometric functions
	xplode$define.pf = define.pf
  
  #(the names of the parameters in the design matrix)
  #to do: change from factor to treatment
	if(!is.na(name.factor)){
	  xplode$factor.col = which(names(model.frame(model)) == name.factor)
	  xplode$factor.colname = name.factor
	  n.factors = length(levels(model.frame(model)[,xplode$factor.col]))
	  factor.parnames = vector("character", n.factors)
	  for(i in 1:n.factors){
	    xplode$factor.parnames[i] = paste(name.factor, levels(model.frame(model)[,xplode$factor.col])[i], sep = "")
	  }
	}
  #name of the continuous predictor
	xplode$cont.col = which(names(model.frame(model)) == name.cont)
	if(!is.na(xplode$cont.col)){
	  xplode$cont.colname = name.cont
	}else{
	  print("Warning: The name of the continuous predictor is not among the names of the model.frame")
	}
  
  
	#re-arrange the fixed effects into n psychometric functions
	#rename the n psychometric functions
	xplode$psychometrics = vector(mode = "list", length = xplode$n.psych.fun)
	names.pf = character(length = xplode$n.psych.fun)
	
  for(i in 1:xplode$n.psych.fun){
		names.pf[i] =  paste("pf", i, sep = "")

    intercept.pointers = define.pf[[i]]$intercept
 		slope.pointers = define.pf[[i]]$slope
    
		#estimate and Variance of the Intercept
		xplode$psychometrics[[i]]$intercept = numeric(length = 2)
    names(xplode$psychometrics[[i]]$intercept) = c("Estimate", "Variance")
    xplode$psychometrics[[i]]$intercept[1] = sum(fixef(model)[intercept.pointers])
    
    if(length(intercept.pointers) == 1)
      {
      xplode$psychometrics[[i]]$intercept[2] = xplode$fixef.vcov[intercept.pointers, intercept.pointers]
    }else
      {
      kombo.intercept = kombo(intercept.pointers)
      xplode$psychometrics[[i]]$intercept[2] =
        xplode$fixef.vcov[kombo.intercept[[1]][1,1], kombo.intercept[[1]][2,1]] +
        xplode$fixef.vcov[kombo.intercept[[1]][1,2], kombo.intercept[[1]][2,2]] +
        (2 * xplode$fixef.vcov[kombo.intercept[[2]][1,1], kombo.intercept[[2]][2,1]])
    }
		
		#estimate and Variance of the Slope
		xplode$psychometrics[[i]]$slope = numeric(length = 2)
		names(xplode$psychometrics[[i]]$slope) = c("Estimate", "Variance")
		xplode$psychometrics[[i]]$slope[1] = sum(fixef(model)[slope.pointers])
		if(length(slope.pointers) == 1){
		  xplode$psychometrics[[i]]$slope[2] = xplode$fixef.vcov[slope.pointers, slope.pointers]
		}else{
		  kombo.slope = kombo(slope.pointers)
		  xplode$psychometrics[[i]]$slope[2] =
		    xplode$fixef.vcov[kombo.slope[[1]][1,1], kombo.slope[[1]][2,1]] +
		    xplode$fixef.vcov[kombo.slope[[1]][1,2], kombo.slope[[1]][2,2]] +
		    (2 * xplode$fixef.vcov[kombo.slope[[2]][1,1], kombo.slope[[2]][2,1]])
		}
	}
	names(xplode$psychometrics) = names.pf
  xplode$psychometrics$pf1$cov = xplode$fixef.vcov[name.cont,"(Intercept)"]
	
  #2)random effects (univariate or multivariate)-------------
	#Str. Dev.
	xplode$ranef.stddev <- as.numeric(attr(VarCorr(model)[[1]], "stddev"))
	
	#Variance-Covariance Matrix of the random effects (two random effect and correlation)
	if(length(xplode$ranef.stddev) > 1){   
		xplode$ranef.VarCov = nearPD(VarCorr(model)[[1]])$mat
		xplode$multirand = TRUE	
		}
				
	#3)family
	#xplode$family = family(model)
  #Bug fix 05.09.2014:
	xplode$family$family = summary(model)$family
	xplode$family$link = summary(model)$link
  
	#number of subjects and names
	xplode$Gp = getME(model, name = "Gp")
	xplode$Groups.levels = levels(getME(model, name = "flist")[[1]])
	xplode$flist = getME(model, name = "flist")
  
	#find the column number in the frame specifiyng for the grouping factor
	xplode$Groups.colnames = names(xplode$flist)
	
  xplode$Groups = length(xplode$Groups.levels)

  xplode$nobs = nobs(model)
	xplode$size = rowSums(model.response(model.frame(model))) #----------------------------> check this with binaries 

  #03.01.14: Check this
  if(is.na(names.response)){
    xplode$response.colnames = colnames(model.response(model.frame(model))) 
	}else{
	  xplode$response.colnames = names.response
	}

  xplode$model.frame = model.frame(model)
  xplode$model.matrix = model.matrix(model)
  
	#Rows per each Subject
	xplode$rps = table(getME(model, "flist")[[1]])
  
  #formula
  xplode$formula = formula(model)
    
  #modified on 14.11.2013 (yet not tested)	
	class(xplode) = "xplode"
	
	#create the methds for the class, i.e., print, summary, ect.
	#print.xplode <- function(x, ...){
	#		cat("This is my vector:\n")
	#		cat(paste(x[1:5]), "...\n")
	#}
	
	return(xplode)
  }