predict_sens <- function(sppp_code, 
                         int_int, int_cwd, int_cwd2, int_pet, int_pet2, 
                         cwd_int, cwd_cwd, cwd_cwd2, cwd_pet, cwd_pet2, 
                         pet_int, pet_cwd, pet_cwd2, pet_pet, pet_pet2){
  ## Function used to predict species' sensitivity rasters based on historic 
  ## climate and randomly drawn parameters from second stage model.
  
  select <- dplyr::select
  
  sp_df <- (sp_clim %>% 
              filter(sp_code == sppp_code) %>% 
              pull(clim_historic_sp))[[1]] %>% 
    lazy_dt()
  
  sp_df <- sp_df %>% 
    rename(cwd_hist = cwd,
           pet_hist = pet) %>%
    mutate(cwd_sens = cwd_int + (cwd_cwd * cwd_hist) + (cwd_cwd2 * cwd_hist * cwd_hist) + (cwd_pet * pet_hist) + (cwd_pet2 * pet_hist * pet_hist),
           pet_sens = pet_int + (pet_cwd * cwd_hist) + (pet_cwd2 * cwd_hist * cwd_hist) + (pet_pet * pet_hist) + (pet_pet2 * pet_hist * pet_hist),
           intercept = int_int + (int_cwd * cwd_hist) + (int_cwd2 * cwd_hist * cwd_hist) + (int_pet * pet_hist) + (int_pet2 * pet_hist * pet_hist)) %>% 
    select(-cwd_hist, -pet_hist) %>% 
    as_tibble()
  
  return(sp_df)
}