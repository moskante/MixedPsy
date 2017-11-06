#' PSE/JND for Multivariable GLMM Using Delta Methods
#'
#' Estimate the Point of Subjective Equivalence (PSE), the Just Noticeable
#' Difference (JND) and the related Standard Errors for a multivariate distribution by means of Delta Method.
#' The method applies to multivariable GLMM having a \emph{probit} link function.
#' The function is based on a recursive use of \code{glmer} and
#' \code{MixDelta}
#'
#' @param xplode.obj an object of class \code{xplode.obj}. The fitted model
#' (object of class \code{"\linkS4class{merMod}"}) from \code{xplode.obj} includes
#' one continuous predictor and one factorial predictor.
#' @param datafr  the data frame fitted with the GLMM model
#'
#' @details The function \code{MixTreatment} is based on a recursive use of
#' \code{glmer} and \code{PsychDelta} to multivariable GLMM including
#' continuous and factorial predictors. The same caveats of \code{PsychDelta}
#' apply (e.g., confidence interval based on normality assumption).
#'
#' @return A list, whose lenght is equal to the levels of the factorial predictor, i.
#' Each cell of the list is equal to the output of \code{delta.psy.probit} applied to
#' a multivariable model whose baseline is level i of the factorial predictor.
#'
#' @references
#' Moscatelli, A., Mezzetti, M., & Lacquaniti, F. (2012). Modeling psychophysical data 
#' at the population-level: The generalized linear mixed model. 
#' Journal of Vision, 12(11):26, 1-17. https://doi.org/10.1167/12.11.26
#'
#' @seealso \code{\link[lme4]{glmer}} for Generalized Linear Mixed Models (including
#' random effects).\code{\link{MixDelta}} for univariable model with delta method.
#' \code{\link{pseMer}} for bootstrap-based confidence intervals.
#'
#' @keywords DeltaMethod Multivariable GLMM
#'
#' @examples
#' library(lme4)
#' data(vibro_exp3)
#' formula.mod <- cbind(faster, slower) ~ speed * vibration + (1 + speed| subject)
#' mod <- glmer(formula = formula.mod, family = binomial(link = "probit"), data = vibro_exp3)
#' xplode.mod <- xplode(model = mod, name.cont = "speed", name.factor = "vibration")
#' MixTreatment(xplode.mod, vibro_exp3)
#'
#' @importFrom stats binomial contrasts<- contr.treatment
#' @export
#'
MixTreatment <- function(xplode.obj, datafr) {

    treat.lev = nlevels(xplode.obj$model.frame[, xplode.obj$factor.col])
    temp.models = delta.par = temp.xplode = vector("list", treat.lev)
    names(delta.par) = xplode.obj$factor.parnames

    for (i in 1:treat.lev) {
        contrasts(datafr[, which(names(datafr) == xplode.obj$factor.colname)]) = contr.treatment(treat.lev,
            base = i)
        temp.models[[i]] = glmer(formula = xplode.obj$formula, family = binomial("probit"), data = datafr,
            nAGQ = 1)
        temp.xplode[[i]] = xplode(temp.models[[i]], name.cont = xplode.obj$cont.colname, name.factor = xplode.obj$factor.colname,
            define.pf = xplode.obj$define.pf)
        delta.par[[i]] = MixDelta(temp.xplode[[i]])[[1]]
    }
    return(delta.par)
}
