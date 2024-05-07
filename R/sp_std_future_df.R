sp_std_future_df <- function(cmip_df, hist_clim_vals, pet_mean, pet_sd, cwd_mean, cwd_sd, temp_mean, temp_sd){
  valid_locations <- hist_clim_vals %>% select(x,y)
  cmip_df <- valid_locations %>% 
    left_join(cmip_df, by = c("x", "y"))
  cmip_df <- cmip_df %>% 
    mutate_at(vars(starts_with("cwd")), 
              ~sp_standardize(.x, cwd_mean, cwd_sd)) %>% 
    mutate_at(vars(starts_with("pet")), 
              ~sp_standardize(.x, pet_mean, pet_sd)) %>% 
    mutate_at(vars(starts_with("temp")), 
              ~sp_standardize(.x, temp_mean, temp_sd))
  return(cmip_df)
}