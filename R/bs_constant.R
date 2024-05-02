bs_constant <- function(data){
  # data <- data %>% # Needed to add this since block bootstrap is returning nested tibble
  #   unnest(cols = c(data))
  const_sens <- data %>% 
    summarise(cwd_const_sens = weighted.mean(cwd_coef, cwd_errorweights),
              pet_const_sens = weighted.mean(pet_coef, pet_errorweights),
              int_const_sens = weighted.mean(int_coef, pet_errorweights))
  return(const_sens)
}