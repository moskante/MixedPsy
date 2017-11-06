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
#' @seealso PsySimulate() for simulating dataframes.
#'
#' @usage data(simul_data)
"simul_data"
#
#' Data from tactile discrimination task - EXP3
#'
#' A dataset containing the response and stimuli from a tactile discrimination task (nine participants). In a
#'  forced-choice experiment, participants were required to discriminate the motion speed of a moving surface by
#'  touching it. Simultaneously with the motion stimulus, a 32Hz masking vibration occurred in half of the trials. 
#'
#' @format A data frame with 72 rows and 6 variables:
#' \describe{
#' \item{speed}{a numeric vector giving the motion speed in cm/s of the moving surface}
#' \item{vibration}{a numeric vector giving the vibration frequency in Hz of the masking stimulus. Either 32Hz or
#'  0 (no vibration - control condition)}
#' \item{faster}{The proportion of trials in which the comparison stimulus was judged as faster than the reference}
#' \item{slower}{The proportion of trials in which the comparison stimulus was judged as slower than the reference}
#' \item{subject}{subject id}
#' }
#' 
#' @source Original data were published in Dallmann et al. (2015).
#' @references 
#' Dallmann, C. J., Ernst, M. O., & Moscatelli, A. (2015). 
#' The role of vibration in tactile speed perception. Journal of Neurophysiology, 114(6),
#'  3131–3139. <doi:10.1152/jn.00621.2015>
#'
#' @usage data(vibro_exp3)
"vibro_exp3"
