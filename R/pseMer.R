#' PSE/JND for GLMM Using Bootstrap Methods
#'
#' Estimates the Point of Subjective Equivalence (PSE), the Just Noticeable
#' Difference (JND) and the related Standard Errors by means of Bootstrap Method.
#' 
#' @param mer.obj An object of class \code{"\linkS4class{merMod}"}.
#' @param B The number of bootstrap samples.
#' @param FUN An optional, custom made function to specify the required parameters to be estimated.
#' if NULL, \code{pseMer()} will estimate the PSE and the JND of a univariable GLMM.
#' @param alpha Significance level of the confidence interval.
#' @param ci.type A vector of character strings representing the type of intervals required. The value 
#' should be any subset of the values c("norm","basic", "stud", "perc", "bca") or simply "all" which will 
#' compute all five types of intervals. "perc" should be always included for the summary table.
#' @param beep Logical. If TRUE, a "ping" sound alerts that the simulation is complete.
#'
#' @return \code{pseMer} returns a list of length 3 including a summary table (Estimate, Standard Error,
#' Inferior and Superior Confidence Interval of the parameters) and the output of  \code{lme4::bootMer()} 
#' and \code{boot::boot.ci()} functions, for further analises. Confidence Intervals in the summary table are
#' based on the percentile method.
#'
#' @details \code{pseMer} estimates PSE and JND (and additional user defined paremters) from a 
#' fitted GLMM model (class \code{"\linkS4class{merMod}"}).
#' 
#' @note A first custom function was written in 2012 for the non-CRAN package MERpsychophisics,
#' based on the algorithm in Moscatelli et al, (2012). The current function is a simple wrapper
#' of \code{lme4::bootMer()} and \code{boot::boot.ci()} functions.
#'
#' @references
#' Moscatelli A, Mezzetti M, Lacquaniti F (2012). Modelling Psychophysical Data
#' at the Population-Level: The Generalized Linear Mixed Model.
#' Journal of Vision, 12(11):26, 1-17.
#' 
#' Bates, D., Maechler, M., Bolker, B., & Walker, S. (2015). Fitting Linear Mixed-Effects
#'  Models Using lme4. Journal of Statistical Software, 67(1), 1–51.
#'
#' @seealso
#' \code{\link{bootMer}} from lme4 package and \code{\link{boot.ci}} from boot package. The "ping" sound is 
#' provided by beepr package.
#'
#' @examples
#' #load simulated data
#' data(psych)
#' formula.mod <- cbind(Longer, Total - Longer) ~ X + (1 + X| Subject)
#' mod1 <- glmer(formula = formula.mod, family = binomial(link = "probit"), data = psych)
#' pse.boot <- pseMer(mod1, B = 600, ci.type = c("norm", "perc"))
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
