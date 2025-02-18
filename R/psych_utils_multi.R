#' Fit Multiple Psychometric Functions with Generalized Linear Models (GLM)
#' 
#' The function fits psychometric functions to data using \code{\link[stats]{glm}} for multiple groups. It supports the use of a binomial error distribution. 
#'  
#' @param data a data frame containing the variables to be used in the model. 
#' @param group_factors a character vector specifying the grouping variables in the dataset. If \code{NULL}, the model will be fit to the entire dataset without grouping.
#' @param formula the formula of the \code{\link[stats]{glm}} model. The response should consist of a binomial outcome (e.g., \code{cbind(yes, no)}).
#' @param link the link function. A character string specifying the link function to be used. By default, \code{"probit"} is used. See \code{\link[stats]{glm}} for available link functions. 
#' 
#' @details
#' This function allows the fitting of psychometric functions to grouped data. If grouping variables are provided through \code{group_factors}, separate models are fit to each group. The function returns a list of models, one for each group, where the model for each group is fitted using the specified \code{formula} and \code{link}.
#' 
#' The models are returned as a named list, with each list element containing the fitted GLM model and the associated group-level information.
#' 
#' 
#' @seealso \code{\link[stats]{glm}}, \code{\link{PsychParameters}}
#' 
#' @examples
#' model_list <- PsychModels(formula = cbind(Longer, Total - Longer) ~ X,
#' data = simul_data,
#' group_factors = "Subject")
#' 
#' model_list_vibro <- PsychModels(vibro_exp3,
#' group_factors = c("subject", "vibration"),
#' formula = cbind(faster, slower) ~ speed)
#' 
#' @importFrom dplyr group_by group_split across distinct
#' @importFrom purrr map_chr map
#' @importFrom magrittr "%>%"
#' @importFrom rlang set_names
#' @importFrom tidyselect all_of
#' 
#' @export 
#'
#' 

PsychModels <- function(data, group_factors = NULL, formula, link = "probit") {
  
  if(!is.null(group_factors)){
    
    grouped_data <- data %>% 
      group_by(across(all_of(group_factors))) %>%
      group_split() 
    group_names <- map_chr(grouped_data, function(group_data) {
      paste(unlist(distinct(group_data, across(all_of(group_factors))) %>% as.list()), collapse = "_")
    })
    
    models <- grouped_data %>%
      set_names(group_names) %>%
      map(function(group_data) {
        group_levels <- distinct(group_data, across(all_of(group_factors))) %>% as.list()
        model <- glm(formula = formula, 
                     family = binomial(link = link), 
                     data = group_data)
        list(model = model, group = group_levels)
      })
    
  } else {
    model <- glm(formula = formula, 
                 family = binomial(link = link), data = data)
    models[["all"]] <- list(model = model, group = "all")
  }
  
  return(models)
}

#' Calculate PSE and JND Parameters from a List of GLM Models
#' 
#' This function calculates the Point of Subjective Equality (PSE) and Just Noticeable Difference (JND) from a list of fitted Generalized Linear Models (GLMs). It extracts these parameters using the \code{\link{PsychDelta}} function and returns them in a structured dataframe.
#' 
#' @param model_list A structured list of grouped models obtained from \code{\link{PsychModels}}. The function can also take as input a GLM model or a list of GLM models.  
#' @param se Logical. if \code{TRUE}, the function includes columns for standard errors of JND and PSE. Default is \code{TRUE}.
#' 
#' @details
#' The function supports three types of input:
#' \itemize{
#'   \item A structured list of models (as produced by \code{\link{PsychModels}}): Extracts PSE and JND for each model and includes the corresponding grouping factors in the output.
#'   \item A single GLM model: Returns a one-row data frame with PSE, JND, and (if requested) standard errors.
#'   \item A list of GLM models: Computes PSE and JND for each model and returns a data frame.
#' }
#' 
#' @return A data frame containing PSE and JND estimates, along with their standard errors (if \code{se = TRUE}). 
#' If the input is a grouped list of models, the output includes columns for the grouping factors.
#' 
#' @seealso \code{\link{PsychModels}}, \code{\link{PsychDelta}}
#' 
#' @examples
#' model_list <- PsychModels(formula = cbind(Longer, Total - Longer) ~ X,
#' data = simul_data,
#' group_factors = "Subject")
#' psych_parameters <- PsychParameters(model_list)
#' 
#' model_list_vibro <- PsychModels(vibro_exp3,
#' group_factors = c("subject", "vibration"),
#' formula = cbind(faster, slower) ~ speed)
#' psych_parameters_vibro <- PsychParameters(model_list_vibro)
#'  
#' @importFrom purrr map_dfr
#' @export
#'
#' 

PsychParameters <- function(model_list, se = TRUE) {
  
  # Helper function to extract parameters
  extract_parameters <- function(model) {
    psejnd <- PsychDelta(model)
    
    # Extract required parameters
    params <- list(
      pse = psejnd["pse", "Estimate"],
      jnd = psejnd["jnd", "Estimate"]
    )
    if (isTRUE(se)) {
      params$pse_se <- psejnd["pse", "Std. Error"]
      params$jnd_se <- psejnd["jnd", "Std. Error"]
    }
    return(params)
  }
  
  
  # Process single model
  if (inherits(model_list, "glm")) {
    return(as.data.frame(extract_parameters(model_list)))
  }
  
  # Process multiple models
  if(all(sapply(model_list, function(x) inherits(x, "glm")))){
    parameters <- map(model_list, extract_parameters) %>%
      dplyr::bind_rows() 
    return(parameters)
  }
  
  if (is.list(model_list) && all(sapply(model_list, function(x) is.list(x) && "model" %in% names(x) && "group" %in% names(x)))) {
    
    parameters <- map_dfr(names(model_list), function(group_key) {
      model_entry <- model_list[[group_key]]  
      params <- extract_parameters(model_entry$model)  
      params_df <- as.data.frame(params)  
      
      # Add group factor columns
      group_info <- model_entry$group
      for (factor_name in names(group_info)) {
        params_df[[factor_name]] <- group_info[[factor_name]]
      }
      
      return(params_df)
    })
    
    return(parameters)
  }
  
  stop("Invalid input: Expected a glm, glmerMod, or a list of such models from PsychModels.")
  
}

#' Interpolate Predictions from a List of GLM Models
#' 
#' This function generates an interpolated dataset by predicting values across a range of an independent variable from a list of generalized linear models (GLMs).
#' 
#' @param model_list A structured list of grouped models obtained from \code{\link{PsychModels}}.  
#' @param n_points An integer number. It specifies the number of points to interpolate along the independent variable range. Default is 100.
#' 
#' @details
#' The function takes a structured list of models, as produced by \code{\link{PsychModels}}, and generates a new dataset with interpolated values for the independent variable. 
#' Predictions are computed at evenly spaced points across the observed range for each model, and the results are returned in a long-format data frame.
#' 
#' @return A data frame containing the interpolated independent variable, the corresponding predicted values from the GLM model, and columns for the grouping factors.
#' 
#' @seealso \code{\link{PsychModels}}, \code{\link[stats]{predict}}.
#' 
#' @examples
#' model_list <- PsychModels(formula = cbind(Longer, Total - Longer) ~ X,
#' data = simul_data,
#' group_factors = "Subject")
#' 
#' longData <- PsychInterpolate(model_list)
#' 
#' # use the interpolated dataset to plot model:
#' library(ggplot2)
#' ggplot(longData, aes(X, prediction, color = Subject)) +
#' geom_line() +
#' geom_point(data = simul_data, aes(X, Longer/Total))
#' 
#' @export
#'
#' 

PsychInterpolate <- function(model_list, n_points = 100) {
  
  if (is.list(model_list) && all(sapply(model_list, function(x) is.list(x) && "model" %in% names(x) && "group" %in% names(x)))){
    # Extract group factors from the first model entry
    first_group_info <- model_list[[1]]$group
    group_factors <- names(first_group_info)  
    
    # Extract predictor variable from the first model formula
    first_model <- model_list[[1]]$model
    formula_terms <- all.vars(formula(first_model))
    x_var <- formula_terms[3]  
    
    # Generate interpolated data for each model in the list
    longData <- map_dfr(names(model_list), function(group_key) {
      model_entry <- model_list[[group_key]]  # Get the model and group info
      model <- model_entry$model
      group_info <- model_entry$group
      
      # Generate evenly spaced values over x range
      x_range <- range(model$model[[x_var]], na.rm = TRUE)
      interp_x <- data.frame(seq(x_range[1], x_range[2], length.out = n_points))
      colnames(interp_x) <- x_var
      
      # Add group factor columns
      for (factor_name in names(group_info)) {
        interp_x[[factor_name]] <- group_info[[factor_name]]
      }
      
      interp_x$prediction <- predict(model, newdata = interp_x, type = "response")
      
      return(interp_x)
    })
    
    return(longData)
  }
  stop("Invalid input: Expected a list of models from PsychModels.")
}


