
#' Internal Function: Set global variable
#' 
#' Setting a global variable is necessary in order to correctly define the function to fit with gnlm. 
#' 
#' @param data dataset
#' @param stimuli A string of the stumuli variable.
#' 
#' @seealso [PsychFunction_gnlm()]
#' 
setGlobalVar <- function(data, stimuli) {
  pkg.globals <- new.env()
  
  pkg.globals$x_values <- data[,stimuli]
}

#' Internal Function: Delete global variable
#' 
#' Remove the global variable
#' 
#' 
rmGlobalVar <- function(){
  rm("x_values", envir = .GlobalEnv)
}