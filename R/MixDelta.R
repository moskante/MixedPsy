#' PSE/JND for Univariable GLMM Using Delta Methods
#'
#' Estimate the Point of Subjective Equivalence (PSE), the Just Noticeable
#' Difference (JND) and the related Standard Errors for an univariate distribution 
#' by means of Delta Method.
#' 
#' @details \code{MixDelta} estimates PSE and JND of a univariable psychometric
#' function (object of class \code{"glm"}).The method only applies to univariable GLMMs 
#'  having a \emph{probit} link function. Use \code{MixTreatment} for multivariable GLMMs.
#'
#' @param xplode.obj an object of class \code{xplode.obj} (univariable GLMMs).
#' @param alpha significance level of the confidence interval. Default is 0.05.
#'
#' @return \code{MixDelta} returns a list of length 1 including Estimate, Standard Error,
#' Inferior and Superior Confidence Interval of PSE and JND. Confidence Intervals
#' are computed as: \eqn{Estimate +/- z(1-(\alpha/2)) * Std.Error}.
#'
#' @note The function assumes that the first model coefficient is the intercept
#' and the second is the slope. The estimate of the JND assumes a \emph{probit}
#' link function.
#'
#' @references
#' Moscatelli, A., Mezzetti, M., & Lacquaniti, F. (2012). Modeling psychophysical data 
#' at the population-level: The generalized linear mixed model. 
#' Journal of Vision, 12(11):26, 1-17. https://doi.org/10.1167/12.11.26
#' 
#' Casella, G., & Berger, R. L. (2002). Statistical inference (2nd ed.). 
#' Pacific Grove, CA: Duxbury Press
#'
#' @seealso
#'  \code{\link{MixTreatment}} for univarible and multivariable GLMM. 
#'  \code{\link{pseMer}} for bootstrap-based confidence intervals. 
#'  \code{\link{xplode}} objects of class \code{xplode.obj}. 
#'
#' @keywords DeltaMethod Univariable GLMM
#' 
#' @examples
#' library(lme4)
#' data(vibro_exp3)
#' formula.mod <- cbind(faster, slower) ~ speed + (1 + speed| subject)
#' mod <- glmer(formula = formula.mod, family = binomial(link = "probit"), 
#'               data = vibro_exp3[vibro_exp3$vibration == 0,])
#' define.mod <- list(pf = list(intercept = 1, slope = 2))
#' xplode.mod <- xplode(model = mod, name.cont = "speed", define.pf = define.mod)
#' pse.jnd <- MixDelta(xplode.mod)
#' 
#' @importFrom stats qnorm
#' @importFrom grDevices palette
#' @export
#'
MixDelta <- function(xplode.obj, alpha = 0.05) {

    # check if link = probit
    if (xplode.obj$family$link != "probit") {
        output = NA
        print("Use a probit link function")
    } else {
        n.pf = length(xplode.obj$psychometrics)
        output = vector("list", length = n.pf)
        names(output) = names(xplode.obj$psychometrics)

        for (i in 1:n.pf) {
            # copy all the variables in temporary objects
            pse <- -(xplode.obj$psychometrics[[i]]$intercept[1]/xplode.obj$psychometrics[[i]]$slope[1])
            slope <- xplode.obj$psychometrics[[i]]$slope[1]

            var.intercept <- xplode.obj$psychometrics[[i]]$intercept[2]
            var.slope <- xplode.obj$psychometrics[[i]]$slope[2]

            # cov(alpha, slope): for all pfs, is approximated to the cov(alpha1, slope1)
            cov.intercept.slope <- xplode.obj$psychometrics$pf1$cov

            # compute all the other variables
            var.pse <- (1/slope^2) * (var.intercept + (2 * pse * cov.intercept.slope) + (pse^2 *
                var.slope))  #PSE
            inferior.pse <- pse - (qnorm(1 - (alpha/2)) * sqrt(var.pse))
            superior.pse <- pse + (qnorm(1 - (alpha/2)) * sqrt(var.pse))

            jnd <- qnorm(0.75) * (1/slope)
            var.jnd <- (qnorm(0.75) * (-1/slope^2))^2 * var.slope  #JND
            inferior.jnd <- jnd - (qnorm(1 - (alpha/2)) * sqrt(var.jnd))
            superior.jnd <- jnd + (qnorm(1 - (alpha/2)) * sqrt(var.jnd))

            output[[i]] <- matrix(rbind(c(pse, sqrt(var.pse), inferior.pse, superior.pse), c(jnd,
                sqrt(var.jnd), inferior.jnd, superior.jnd)), nrow = 2, dimnames = list(param <- c("pse",
                "jnd"), statistics = c("Estimate", "Std. Error", "Inferior", "Superior")))
        }
    }


    return(output)
}
