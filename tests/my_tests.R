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


