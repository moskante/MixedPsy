source("MERpsychophysics.r")

#this ensure that your output will match mine
set.seed(32787)																	                                 

#psychometric function to fit single-subject data
datafr.S1 <- MERsimulate(fixeff = c(-7.2, 0.0875), nsubject = 1, constant = T)     
plot(Longer/Total ~ X, data = datafr.S1)
fit.S1 = psych.function(ps.formula = cbind(Longer, Total - Longer) ~ X, ps.link = "probit",
                        ps.data = datafr.S1, x.range = c(40, 120), ps.lines = T, ps.col = "black")
fit.S1

#this simulates a second condition or a second participant
datafr.S2 <- MERsimulate(fixeff = c(-6.8, 0.0875), nsubject = 1, constant = T)     
points(Longer/Total ~ X, data = datafr.S2, col = "red")
fit.S2 = psych.function(ps.formula = cbind(Longer, Total - Longer) ~ X, ps.link = "probit",
                        ps.data = datafr.S2, x.range = c(40, 120), ps.lines = T, ps.col = "red")
fit.S2

#Analysis of the clustered data: Mixed Models
#In two steps: Simulate a dataset with a categorical variable ("condition")
datafr.1 <- MERsimulate(fixeff = c(-7.5, 0.0875), nsubject = 6, constant = T)                               
levels(datafr.1$Subject) = c("S1", "S2", "S3", "S4", "S5", "S6")
datafr.1$condition = rep("A", 54)

datafr.2 <- MERsimulate(fixeff = c(-7, 0.0875),nsubject = 6, constant = T)                              
levels(datafr.2$Subject) = c("S1", "S2", "S3", "S4", "S5", "S6")
datafr.2$condition = rep("B", 54)

datafr = merge(datafr.1, datafr.2, all = T)
datafr$condition = as.factor(datafr$condition)

#How to estimate the PSE in the two condtions?
#1)use MERtreatment (note, you have to provide also the dataframe to the function)
formula.mod = cbind(Longer, Total - Longer) ~ X * condition + (1 + X| Subject)
mod1 <- glmer(formula = formula.mod, family = binomial(link = "probit"), data = datafr)
xplode.mod1 = xplode.mer(model = mod1, name.cont = "X", name.factor = "condition")
MERtreatment(xplode.mod1, datafr)

#MERtreatment(.) assumes by default that the continuous predictor is the first in the formula. If you use
#a different order, you have to specify which parameter pertains to the intercept and which to the slope
#using define.pf. See the example below, I fit the same model inverting the oredr of predictors
formula.mod = cbind(Longer, Total - Longer) ~ condition * X + (1 + X| Subject)
mod1 <- glmer(formula = formula.mod, family = binomial(link = "probit"), data = datafr)
define.mod = list(pf1 = list(intercept = 1, slope = 3),
                  pf2 = list(intercept = c(1,2), slope = c(3,4)))
xplode.mod1 = xplode.mer(model = mod1, name.cont = "X", name.factor = "condition",
                         define.pf = define.mod)
MERtreatment(xplode.mod1, datafr)

#write your function to estimate jnd and pse with a bootstrap method:
fun2mod1 = function(mer.obj){
  #allocate space: 4 parameters (jnd_A, jnd_B, pse_A, pse_B)
  jndpse = vector(mode = "numeric", length = 4)
  names(jndpse) = c("jnd_A","jnd_B", "pse_A", "pse_B")
  
  jndpse[1] = qnorm(0.75)/fixef(mer.obj)[2] #jnd_A
  jndpse[2] = qnorm(0.75)/(fixef(mer.obj)[2] + fixef(mer.obj)[4]) #jnd_B
  jndpse[3] = -fixef(mer.obj)[1]/fixef(mer.obj)[2] #pse_A
  jndpse[4] = -(fixef(mer.obj)[1] + fixef(mer.obj)[3])/(fixef(mer.obj)[2] + fixef(mer.obj)[4]) #pse_B
  return(jndpse)
}

#Now use the pseMer function
BootEstim = pseMer(mod1, B = 200, FUN = fun2mod1)

#...and the estimates:
BootEstim[[1]]

#quick plotting
library(ggplot2)
datafr$predict.mod1 = predict(mod1, type = "response")
ggplot(datafr, aes(x = X, y = Longer/Total, color = condition)) +
  ylab("Probability") + xlab("X") +
  geom_point() +
  geom_line(aes(x = X, y = predict.mod1)) +
  facet_wrap(~ Subject, ncol = 3)

#2)alternatively, run twice the model with a different baseline
#This way we can choose the baseline 
datafr$conditionA = ifelse(datafr$condition == "A", 0 , 1)

#fit the GLMM(probit link function)
formula.modA = cbind(Longer, Total - Longer) ~ X * conditionA + (1 + X| Subject)
modA <- glmer(formula = formula.modA, family = binomial(link = "probit"), data = datafr, nAGQ = 1)
summary(modA)

#New: xplode.mer (indicate the names of the predictor)
xplode.modA = xplode.mer(model = modA, name.cont = "X", name.factor = "conditionA")

#run MERdelta.probit on the xplode object
#This is the PSE and JND for A (B is coded as a dummy == 1)
MERdelta.probit(xplode.modA)

#pse and jnd CI with the bootstrap (this takes time). The function based is a very simple one, based on bootMer{lme4.1}
estimA = pseMer(modA)

#Now condition B is the baseline
datafr$conditionB = ifelse(datafr$condition == "B", 0 , 1)

#fit the GLMM(probit link function)
formula.modB = cbind(Longer, Total - Longer) ~ X * conditionB+ (1 + X| Subject)
modB <- glmer(formula = formula.modB, family = binomial(link = "probit"), data = datafr, nAGQ = 1)

xplode.modB = xplode.mer(model = modB, name.cont = "X", name.factor = "conditionB")
#This is the PSE and JND for B (A is coded as a dummy == 1)
MERdelta.probit(xplode.modB)

#pse and jnd CI with the bootstrap (this takes time). The function based is a very simple one, based on bootMer{lme4.1}
estimB = pseMer(modB)

#3) I would like to implement a third method based on the variance-covariance matrix.
#For the moment is not ready, as I don't know how to find the three-terms covariance. Without that, the CI is unprecise.
#Define the psychometric function and run the model just once
define.mod = list(pf1 = list(intercept = 1, slope = 2),
                    pf2 = list(intercept = c(1,3), slope = c(2,4)))

xplode.mod = xplode.mer(model = modA, name.cont = "X", name.factor = "condition", define.pf = define.mod)
MERdelta.probit(xplode.mod)

