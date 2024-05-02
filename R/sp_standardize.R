sp_standardize <- function(val, sp_mean, sp_sd){
  std_val <- (val - sp_mean) / sp_sd
  return(std_val)
}