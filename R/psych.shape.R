#psych.shape() plot a psychometric function of a knwonwn shape and location.
#The function always requires an existing plot (add = TRUE)
#Arguments
#pse, jnd  the pse and the jnd of the desired psychometric function
#x.range  the range over which the function will be plotted.
#ps.link		character. link function
#ps.col     character. Color of the fitted function
#ps.lwd     line width
#Example
#y = c(0,1)
#x = c(-40, 40)
#plot(y ~ x, type = "n", bty = "n", lab = c(5,3,7))
#psych.shape(pse = 0, jnd = 6, x.range = c(-40, 40), ps.col = "gray", ps.lwd = 3)
#psych.shape(pse = 6, jnd = 6, x.range = c(-40, 40), ps.col = "black")
#psych.shape(pse = 6, jnd = 6, x.range = c(-40, 40), ps.col = "red", ps.link = "logit", ps.lwd = 3)

psych.shape = function(pse = 0, jnd, x.range = c(NA, NA),  ps.link = "probit", ps.col = "black",
                       ps.lwd = 1){
  if(ps.link == "probit"){
    slope = qnorm(0.75) * (1/jnd)
    curve(expr = pnorm(x, mean = pse, sd = 1/slope), from = x.range[1], to = x.range[2], col = ps.col,
          add = TRUE, lwd = ps.lwd)
  }else{
    if(ps.link == "logit"){
      slope = log(3) * (1/jnd)
      curve(expr = plogis(x, location = pse, scale = 1/slope), from = x.range[1], to = x.range[2], col = ps.col,
            add = TRUE, lwd = ps.lwd)
    }else{
    warning("For the moment it only works with probit and logit link function")
      }
    }
  return(0)
  }