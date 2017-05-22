#' A simulated psychophysical dataset
#'
#' A dataset containing simulated data for 8 subjects. Created using 
#' \code{PsySimulate(constant = T)}. The variables are as follows:
#'
#' @format A data frame with 72 rows and 6 variables:
#' \describe{
#' \item{X}{Samples of the X interval \code{c(40,120)}}
#' \item{Intercept}{Intercept of the psychometric function}
#' \item{Slope}{Slope of the psychometric function}
#' \item{Longer}{Number of trials in which response is judged "longer" than standard}
#' \item{Total}{Total number of trials for each sample of X interval}
#' \item{Subject}{Subject code (S1 to S8)}
#' }
#' 
#' @seealso PsySimulate()
#'
#' @usage data(psych)
"psych"