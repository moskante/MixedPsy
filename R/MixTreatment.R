#' PSE/JND for Multivariable GLMM Using Delta Methods
#'
#' Estimates the Point of Subjective Equivalence (PSE), the Just Noticeable
#' Difference (JND) and the related Standard Errors by means of Delta Method.
#' The method apply to multivariable GLMM having a \emph{probit} link function.
#' The function is based on a recursive use of \code{glmer} and
#' \code{delta.psy.probit}
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
#' Moscatelli A, Mezzetti M, Lacquaniti F (2012). Modelling Psychophysical Data at the
#' Population-Level: The Generalized Linear Mixed Model.
#' Journal of Vision, 12(11):26, 1-17.
#'
#' @seealso \code{\link[lme4]{glmer}} for Generalized Linear Mixed Models (including
#' random effects).\code{\link{PsychFunction}} for a univariable model.
#' \code{\link{pseMer}} provides the bootstrap-based confidence intervals.
#'
#' @keywords Delta Method, Multivariable GLMM #####NOT SHOWN
#'
#' @examples
#' #In two steps: Simulate a dataset with a categorical variable ("condition")
#' datafr.1 <- PsySimulate(fixeff = c(-7.5, 0.0875), nsubject = 6, constant = TRUE)
#' levels(datafr.1$Subject) = c("S1", "S2", "S3", "S4", "S5", "S6")
#' datafr.1$condition = rep("A", 54)
#'
#' datafr.2 <- PsySimulate(fixeff = c(-7, 0.0875),nsubject = 6, constant = TRUE)
#' levels(datafr.2$Subject) = c("S1", "S2", "S3", "S4", "S5", "S6")
#' datafr.2$condition = rep("B", 54)
#'
#' datafr = merge(datafr.1, datafr.2, all = TRUE)
#' datafr$condition = as.factor(datafr$condition)
#'
#' #How to estimate the PSE in the two condtions?
#' #1)use MixTreatment (note, you have to provide also the dataframe to the function)
#' library(lme4)
#' formula.mod = cbind(Longer, Total - Longer) ~ X * condition + (1 + X| Subject)
#' mod1 <- glmer(formula = formula.mod, family = binomial(link = "probit"), data = datafr)
#' xplode.mod1 = xplode(model = mod1, name.cont = "X", name.factor = "condition")
#' MixTreatment(xplode.mod1, datafr)
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
