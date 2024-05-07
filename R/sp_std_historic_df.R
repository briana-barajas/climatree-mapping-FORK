sp_std_historic_df <- function(hist_clim_vals, pet_mean, pet_sd, cwd_mean, cwd_sd, temp_mean, temp_sd){
  hist_clim_vals <- hist_clim_vals %>% 
    mutate_at(vars(starts_with("cwd")), 
              ~sp_standardize(.x, cwd_mean, cwd_sd)) %>% 
    mutate_at(vars(starts_with("pet")), 
              ~sp_standardize(.x, pet_mean, pet_sd)) %>% 
    mutate_at(vars(starts_with("temp")), 
              ~sp_standardize(.x, temp_mean, temp_sd))
  return(hist_clim_vals)
}