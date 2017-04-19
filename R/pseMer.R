#CI of the jndpse and jnd based on bootMer{lme4.1} function

pseMer = function(mer.obj, B = 200, FUN = NULL, beep = F){
  
  boot.flag = require("boot")
  if(boot.flag == F){
    print("Package boot is missing. Install the package to use the function.")}
  
  if(is.null(FUN)){
    myfun<-function(mer.obj){
      jndpse = vector(mode = "numeric", length = 2)
      names(jndpse) = c("JND", "PSE")
      jndpse[1] = qnorm(0.75)/fixef(mer.obj)[2]
      jndpse[2] = -fixef(mer.obj)[1]/fixef(mer.obj)[2]
      return(jndpse)}
  }else{
    myfun = match.fun(FUN)
  }

  np = length(myfun(mer.obj))
  parname = names(myfun(mer.obj))
  if(is.null(parname)){print("Warning messages: Parameters have no names")}
  summary = matrix(NA, nrow = np, ncol = 3, dimnames = list(parname, c("Estimate", "Inferior", "Superior")))
  
  #Warning: these two lines are only valid in lme4.1
  boot.samp <- bootMer(mer.obj, myfun, nsim = B)
  summary[,1] = boot.samp$t0
  
  jndpseconf = vector(mode = "list", length = np)
    
  for(i in 1:np){
    jndpseconf[[i]] <- boot.ci(boot.samp,type=c("norm", "basic", "perc"), index = i)
    print(paste(parname[i], " 95% CI:",jndpseconf[[i]]$perc[4],"  ",jndpseconf[[i]]$perc[5]))
    summary[i,2] = jndpseconf[[i]]$perc[4]
    summary[i,3] = jndpseconf[[i]]$perc[5]
  }
  
  if(beep == T){
    beep.flag = require("beepr")
    if(beep.flag == F){
      print("Package beep is missing. Install the package or set beep = F")}
    beep()
  }
  out = list(summary, boot.samp, jndpseconf)
  return(out)
}