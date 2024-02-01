utils::globalVariables(c("x_values"))
x_values <- NULL
#' Internal Function: Set global variable
#' 
#' Setting a global variable is necessary in order to correctly define the function to fit with gnlm. 
#' 
#' @param data dataset
#' @param stimuli A string of the stumuli variable.
#' 
#' @seealso [PsychFunction_gnlm()]
#' @noRd
#' @return Invisible NULL.
#'
#' @keywords internal

setGlobalVar <- function(data, stimuli) {
  #unlockBinding("x_values", globalenv())
  #globalenv()$x_values <- data[, stimuli]
  #assign("x_values", data[, stimuli], envir = .GlobalEnv)
  #x_values <<- data[, stimuli]
  # Lock the binding after modification
  #lockBinding("x_values", globalenv())
  #invisible(NULL)
}

#' Internal Function: Delete global variable
#' 
#' Remove the global variable
#' 
#' 
rmGlobalVar <- function(){
 rm("x_values", envir = .GlobalEnv)
}