## Function to run second stage model for first stage CWD, 
## PET and Intercept terms and organize coefficients
bs_ss <- function(data){
  # data <- data %>% # Needed to add this since block bootstrap is returning nested tibble
  #   unnest(cols = c(data))
  cwd_mod <- data %>% 
    run_ss(outcome = "cwd_coef")
  pet_mod <- data %>% 
    run_ss(outcome = "pet_coef")
  int_mod <- data %>% 
    run_ss(outcome = "int_coef")
  # return(list("int_mod" = c(int_mod),
  #             "cwd_mod" = c(cwd_mod),
  #             "pet_mot" = c(pet_mod)))
  return(list(int_int = int_mod[1],
              int_cwd = int_mod[2],
              int_cwd2 = int_mod[3],
              int_pet = int_mod[4],
              int_pet2 = int_mod[5],
              cwd_int = cwd_mod[1],
              cwd_cwd = cwd_mod[2],
              cwd_cwd2 = cwd_mod[3],
              cwd_pet = cwd_mod[4],
              cwd_pet2 = cwd_mod[5],
              pet_int = pet_mod[1],
              pet_cwd = pet_mod[2],
              pet_cwd2 = pet_mod[3],
              pet_pet = pet_mod[4],
              pet_pet2 = pet_mod[5]))
}
