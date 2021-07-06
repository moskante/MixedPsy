#' A simulated psychophysical dataset
#'
#' A dataset containing simulated responses to a foced-choice task for eight clusters. Created using 
#' \code{PsySimulate()}. The variables are as follows:
#'
#' @format A data frame with 72 rows (9 observaxions x 8 simulated participants) and 6 variables:
#' \describe{
#' \item{X}{samples in the continuous interval (range \code{c(40,120)})}
#' \item{Intercept,Slope}{simulated participant's coefficients (combination of random and fixed effect)}
#' \item{Longer}{number of trials in which response is judged "longer" than standard}
#' \item{Total}{total number of trials for sample in X}
#' \item{Subject}{simulated participant's identification code (S1 to S8)}
#' }
#' 
#' @seealso \code{\link{PsySimulate}} for simulating dataframes with custom parameters.
#'
#' @usage data(simul_data)
"simul_data"
#
#' Data from tactile discrimination task - (Dallmann et al., 2015).
#'
#' A dataset containing the response recorded from a tactile discrimination task (nine participants). In a
#'  forced-choice experiment, participants were required to discriminate the motion speed of a moving surface by
#'  touching it. Simultaneously with the motion stimulus, a 32Hz masking vibration occurred in half of the trials. 
#'
#' @format A data frame with 126 rows (14 observations x 9 participants) and 5 variables:
#' \describe{
#' \item{speed}{numeric, speed of the moving surface (in cm/s, range \code{c(1,16)})}
#' \item{vibration}{factor representing the vibration frequency of the masking stimulus. Two levels: 32 (vibration in the experimental condition, in Hz) 
#' or 0 (no vibration - control condition)}
#' \item{faster}{proportion of trials in which the comparison stimulus was judged as faster than the reference}
#' \item{slower}{proportion of trials in which the comparison stimulus was judged as slower than the reference}
#' \item{subject}{participant's identification code}
#' }
#' 
#' @source Original data were published in Dallmann et al. (2015).
#' @references 
#' Dallmann, C. J., Ernst, M. O., & Moscatelli, A. (2015). 
#' The role of vibration in tactile speed perception. Journal of Neurophysiology, 114(6),
#' 3131–3139. <doi:10.1152/jn.00621.2015>
#'
#' @usage data(vibro_exp3)
"vibro_exp3"
