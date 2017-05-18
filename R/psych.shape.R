psych.shape = function(pse = 0, jnd, x.range = c(NA, NA), ps.link = "probit", ps.col = "black", ps.lwd = 1) {
    if (ps.link == "probit") {
        slope = qnorm(0.75) * (1/jnd)
        curve(expr = pnorm(x, mean = pse, sd = 1/slope), from = x.range[1], to = x.range[2], col = ps.col, 
            add = TRUE, lwd = ps.lwd)
    } else {
        if (ps.link == "logit") {
            slope = log(3) * (1/jnd)
            curve(expr = plogis(x, location = pse, scale = 1/slope), from = x.range[1], to = x.range[2], 
                col = ps.col, add = TRUE, lwd = ps.lwd)
        } else {
            warning("For the moment it only works with probit and logit link function")
        }
    }
    return(0)
}
