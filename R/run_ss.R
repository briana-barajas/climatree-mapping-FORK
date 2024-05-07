## Function defining main second stage model
run_ss <- function(data, outcome = "cwd_coef"){
  if (outcome == "cwd_coef") {
    error_weights = data$cwd_errorweights
  } else if (outcome == "pet_coef") {
    error_weights = data$pet_errorweights
  } else if (outcome == "int_coef") {
    error_weights = data$int_errorweights
  }
  formula <- as.formula(paste(outcome, " ~ cwd.spstd + I(cwd.spstd^2) + pet.spstd + I(pet.spstd^2)"))
  mod <- lm(formula, data=data, weights = error_weights)
  # mod <- lm(formula, data=data)
  coefs <- mod %>% 
    tidy() %>% 
    pull(estimate) 
  return(coefs)
}