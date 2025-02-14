library(MixedPsy)

test_data <- simul_data[simul_data$Subject == "S1",]

# glm model
glm_model <- PsychFunction(formula= cbind( Longer , Total - Longer ) ~ X, response = c("Longer", "Total - Longer"), stimuli = "X", model = "glm", link = "probit", data = test_data)

summary(glm_model$glm)

PsychDelta(glm_model$glm)
PsychPlot(glm_model$glm)

logit_model <- PsychFunction(response = c("Longer", "Total - Longer"), stimuli = "X", model = "glm", link = "logit", data = test_data)

summary(logit_model$glm)
PsychDelta(logit_model$glm)
PsychPlot(logit_model$glm)


# gnlm
my_model <- PsychFunction(formula= cbind( Longer , Total - Longer ) ~ X, response = c("Longer", "Total - Longer"), stimuli = "X", model = "gnlm", link = "probit", data = test_data, guess = FALSE, lapse = FALSE)

my_model1 <- PsychFunction(response = c("Longer", "Total - Longer"), stimuli = "X", model = "gnlm", link = "logit", data = test_data, guess = TRUE, lapse = TRUE)

my_model2 <- PsychFunction(response = c("Longer", "Total - Longer"), stimuli = "X", model = "gnlm", link = "weibull", data = test_data, guess = FALSE, lapse = FALSE)



boot_psych_try <- function(data, indices){
  boot_data <- data[indices,]
  tryCatch(
    {mod <- PsychFunction(response = c("Longer", "Total - Longer"), stimuli = "X", model = "gnlm", link = "probit", data = boot_data, guess = FALSE, lapse = FALSE)
    return(mod$gnlm_coeff)
    },
    error = function(e){
      return(rep(NA, 2))
    }
  )
}

boot_psych <- function(data, indices){
  boot_data <- data[indices,]
  mod <- PsychFunction(response = c("Longer", "Total - Longer"), stimuli = "X", model = "gnlm", link = "probit", data = boot_data, guess = FALSE, lapse = FALSE)
    return(mod$gnlm_coeff)
}




for (i in 1:num_bootstrap_samples) {
  
  # Call boot with try() function to handle errors
  bootstrap_results[i, ] <- try(
    boot(
      data = test_data,
      statistic = boot_psych,
      R = 1,
    )[["t"]],
    silent = TRUE
  )
}

num_bootstrap_samples <- 10
bootstrap_results <- matrix(NA, nrow = num_bootstrap_samples, ncol = 2)

for (i in 1:num_bootstrap_samples) {
  # Generate bootstrap indices
  boot_indices <- sample(nrow(test_data), replace = TRUE)
  
  # Initialize retry counter
  max_retries <- 3
  retry_count <- 0
  
  # Attempt model fitting with retries
  while (retry_count < max_retries) {
    # Call boot with try() function to handle errors
    result <- try(
      boot(
        data = test_data,
        statistic = boot_psych,
        R = 1,
      ),
      silent = TRUE
    )
    
    # Check if model fitting was successful
    if (!inherits(result, "try-error")) {
      # Store the result if successful
      bootstrap_results[i, ] <- result
      break  # Exit retry loop
    } else {
      # Increment retry count
      retry_count <- retry_count + 1
    }
  }
}



# Extract parameter estimates
parameter_estimates <- t(apply(bootstrap_results, 1, function(x) x$t))


