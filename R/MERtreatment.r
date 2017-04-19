MERtreatment = function(xplode.obj, datafr){

  treat.lev = nlevels(xplode.obj$model.frame[,xplode.obj$factor.col])
  temp.models = delta.par = temp.xplode = vector("list", treat.lev)
  names(delta.par) = xplode.obj$factor.parnames
  
  for(i in 1:treat.lev){
   contrasts(datafr[, which(names(datafr) == xplode.obj$factor.colname)]) = contr.treatment(treat.lev, base = i)
   temp.models[[i]] = glmer(formula = xplode.obj$formula, family = binomial("probit"), data = datafr, nAGQ = 1)
   temp.xplode[[i]] = xplode.mer(temp.models[[i]], name.cont = xplode.obj$cont.colname,
                                 name.fact = xplode.obj$factor.colname, define.pf = xplode.obj$define.pf)
   delta.par[[i]] = MERdelta.probit(temp.xplode[[i]])[[1]]
 }
  return(delta.par)
}