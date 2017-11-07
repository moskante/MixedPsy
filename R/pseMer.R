#' PSE/JND for GLMM Using Bootstrap Methods
#'
#' Estimates the Point of Subjective Equivalence (PSE), the Just Noticeable
#' Difference (JND) and the related Standard Errors by means of Bootstrap Method.
#' 
#' @param mer.obj An object of class \code{"\linkS4class{merMod}"}.
#' @param B integer: the number of bootstrap samples.
#' @param FUN An optional, custom made function to specify the required parameters to be estimated.
#' if NULL, \code{pseMer()} will estimate the PSE and the JND of a univariable GLMM.
#' @param alpha Significance level of the confidence interval.
#' @param ci.type A vector of character strings representing the type of intervals required. The value 
#' should be any subset of the values c("norm","basic", "stud", "perc", "bca") or simply "all" which will 
#' compute all five types of intervals. "perc" should be always included for the summary table.
#' @param beep Logical. If TRUE, a "ping" sound alerts that the simulation is complete.
#'
#' @return \code{pseMer} returns a list of length 3 including a summary table (Estimate, Standard Error,
#' Inferior and Superior Confidence Interval of the parameters) and the output of  \code{\link[lme4]{bootMer}} 
#' and \code{\link[boot]{boot.ci}} functions, for further analises. Confidence Intervals in the summary table are
#' based on the percentile method.
#'
#' @details \code{pseMer} estimates PSE and JND (and additional user defined paremters) from a 
#' fitted GLMM model (class \code{"\linkS4class{merMod}"}). 
#' The "ping" sound is provided by \code{\link[beepr]{beep}} function from the \code{beepr} package.
#' 
#' @note A first custom function was written in 2012 for the non-CRAN package MERpsychophisics,
#' based on the algorithm in Moscatelli et al. (2012). The current function is a simple wrapper
#' of \code{lme4::bootMer()} and \code{boot::boot.ci()} functions.
#' 
#' Increasing the nuber of bootstrap samples (\code{B}) makes the estimate more reliable. 
#' However, this will also increase the duration of the computation.
#'
#' @references
#' Moscatelli, A., Mezzetti, M., & Lacquaniti, F. (2012). Modeling psychophysical data 
#' at the population-level: The generalized linear mixed model. 
#' Journal of Vision, 12(11):26, 1-17. https://doi.org/10.1167/12.11.26
#' 
#' Bates, D., Mächler, M., Bolker, B., & Walker, S. (2015). Fitting Linear Mixed-Effects 
#' Models Using lme4. Journal of Statistical Software, 67(1), 51. https://doi.org/10.18637/jss.v067.i01
#'
#' @seealso
#' \code{\link[lme4]{bootMer}} from \code{lme4} package and \code{\link[boot]{boot.ci}} from \code{boot} package. 
#' 
#' @keywords Univariable Multivariable GLMM Bootstrap
#'
#' @examples
#' ## Example 1: estimate pse/jnd of a univariable GLMM
#' library(lme4)
#' data(vibro_exp3)
#' formula.mod1 <- cbind(faster, slower) ~ speed + (1 + speed| subject)
#' mod1 <- glmer(formula = formula.mod1, family = binomial(link = "probit"), 
#'               data = vibro_exp3[vibro_exp3$vibration == 0,])
#' BootEstim.1 <- pseMer(mod1, B = 100, ci.type = c("norm", "perc"))
#' 
#' ## Example 2: specify custom parameters for bootstrap estimation of a 
#' # multivariate model
#' 
#' formula.mod2 <- cbind(faster, slower) ~ speed * vibration + (1 + speed| subject)
#' mod2 <- glmer(formula = formula.mod2, family = binomial(link = "probit"), 
#'                data = vibro_exp3)
#'               
#' fun2mod = function(mer.obj){
#' #allocate space: 4 parameters (jnd_0Hz, jnd_32Hz, pse_0Hz, pse_32Hz) j
#' jndpse = vector(mode = "numeric", length = 4)
#' names(jndpse) = c("jnd_0Hz","jnd_32Hz", "pse_0Hz", "pse_32Hz")
#' jndpse[1] = qnorm(0.75)/fixef(mer.obj)[2] #jnd_0Hz
#' jndpse[2] = qnorm(0.75)/(fixef(mer.obj)[2] + fixef(mer.obj)[4]) #jnd_32Hz
#' jndpse[3] = -fixef(mer.obj)[1]/fixef(mer.obj)[2] #pse_0Hz
#' jndpse[4] = -(fixef(mer.obj)[1] + fixef(mer.obj)[3])/(fixef(mer.obj)[2] 
#'                + fixef(mer.obj)[4]) #pse_32Hz
#' return(jndpse)
#' }
#' 
#' BootEstim.2 = pseMer(mod2, B = 100, FUN = fun2mod)
#' 
#' @export
#' @importFrom lme4 bootMer
#' @importFrom Matrix nearPD
#' @importFrom boot boot.ci
#' @importFrom beepr beep

pseMer <- function(mer.obj, B = 200, FUN = NULL, alpha = 0.05, 
                   ci.type = c("norm", "basic", "perc"), beep = F) {

    if (is.null(FUN)) {
        myfun <- function(mer.obj) {
            jndpse = vector(mode = "numeric", length = 2)
            names(jndpse) = c("JND", "PSE")
            jndpse[1] = qnorm(0.75)/fixef(mer.obj)[2]
            jndpse[2] = -fixef(mer.obj)[1]/fixef(mer.obj)[2]
            return(jndpse)
        }
    } else {
        myfun = match.fun(FUN)
    }

    np = length(myfun(mer.obj))
    parname = names(myfun(mer.obj))
    if (is.null(parname)) {
        print("Warning messages: Parameters have no names")
    }
    summary = matrix(NA, nrow = np, ncol = 3,
                     dimnames = list(parname, c("Estimate", "Inferior","Superior")))

    boot.samp <- lme4::bootMer(mer.obj, myfun, nsim = B)
    summary[, 1] = boot.samp$t0

    jndpseconf = vector(mode = "list", length = np)
    my.conf = 1 - alpha
    
    for (i in 1:np) {
        jndpseconf[[i]] <- boot::boot.ci(boot.samp, conf = my.conf, type = ci.type, index = i)
        print(paste(parname[i], " 95% CI:", jndpseconf[[i]]$percent[4], "  ", jndpseconf[[i]]$percent[5]))
        summary[i, 2] = jndpseconf[[i]]$percent[4]
        summary[i, 3] = jndpseconf[[i]]$percent[5]
    }

    if (beep == T) {
        beepr::beep()
    }
    out = list(summary, boot.samp, jndpseconf)
    return(out)
}
