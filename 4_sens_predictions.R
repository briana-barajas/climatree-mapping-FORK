#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Briana Barajas, Fletcher McConnell, Rosemary Juarez, Vanessa Salgado
# Project: Mapping Tree Species' Drought Sensitivity Under Climate Change
# Institution: Bren School of Environmental Science & Management - UCSB
# Date: 2024-06-07
# Purpose: Predict species growth given estimates for future variations in climatic water deficit through 2100**
#
# Input files:
# - ss_bootstrap (species specific)
# - sp_clim_predictions.gz
# -
#
# Output files (species specific):
# - sp_rwi
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(fixest)
library(sp)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(tidyverse)
library(dtplyr)
library(prediction)
library(tictoc)
library(furrr)
library(snow)
library(profvis)
library(tmap)
library(tidylog)


# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Load data --------------------------------------------------------
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Define path
data_dir <- "~/../../capstone/climatree/raw_data/"
output_dir <- "~/../../capstone/climatree/output/"

for(species in c(spp_code_list)){
  
  n_mc <- 1000
  n_cores <- availableCores()
  future::plan(multisession, workers = n_cores)
  
  my_seed <- 5597
  

# 1. Second stage model
mod_df <- read_rds(paste0(output_dir, "ss_bootstrap_", species, ".rds"))
mod_df <- mod_df %>% 
  rename(iter_idx = boot_id)

# 2. Species-standardized historic and future climate
sp_clim <- read_rds(paste0(output_dir, "sp_clim_predictions.gz")) %>% 
  filter(sp_code == species)
species_list <- sp_clim %>% select(sp_code)
            
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Assign MC coefs and CMIP models  ---------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Join second stage coefficients to species list
sp_mc <- species_list %>% 
  select(sp_code) %>% 
  crossing(iter_idx = seq(n_mc)) %>% 
  left_join(mod_df, by = "iter_idx")


## Assign specific cmip realization to each MC iteration 
n_cmip_mods <- 47
cmip_assignments <- tibble(iter_idx = seq(1, n_mc)) %>% 
  mutate(cmip_idx = sample(seq(n_cmip_mods), n_mc, replace = TRUE))

## Join cmip model assignments
sp_mc <- sp_mc %>% 
  left_join(cmip_assignments, by = "iter_idx")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define functions ---------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


calc_rwi_partials <- function(sppp_code, cmip_id, sensitivity, cwd_const_sens, pet_const_sens, int_const_sens){
  select <- dplyr::select
  
  ## Function used to predict species' RWI rasters based on predicted 
  ## sensitivity raster and assigned CMIP model of future climate. Also
  ## integrates calculations of partialling our climate / sensitivity variations
  
  ## Predict RWI under CMIP scenario
  sp_fut_clim <- sp_clim %>% 
    filter(sp_code == sppp_code) %>% 
    pull(clim_cmip_sp)
  sp_fut_clim <- sp_fut_clim[[1]] %>% 
    lazy_dt() %>% 
    select(x,y,
           cwd_cmip_end = paste0("cwd_cmip_end", as.character(cmip_id)),
           pet_cmip_end = paste0("pet_cmip_end", as.character(cmip_id)),
           cwd_cmip_start = paste0("cwd_cmip_start", as.character(cmip_id)),
           pet_cmip_start = paste0("pet_cmip_start", as.character(cmip_id))) %>% 
    left_join(sensitivity, by = c("x", "y")) %>% 
    mutate(rwi_pred_end = intercept + (pet_cmip_end * pet_sens) + (cwd_cmip_end * cwd_sens),
           rwi_pred_start = intercept + (pet_cmip_start * pet_sens) + (cwd_cmip_start * cwd_sens),
           rwi_pclim_end = int_const_sens + (pet_cmip_end * pet_const_sens) + (cwd_cmip_end * cwd_const_sens),
           rwi_pclim_start = int_const_sens + (pet_cmip_start * pet_const_sens) + (cwd_cmip_start * cwd_const_sens)) %>% 
    select(x,y,
           cwd_sens,
           pet_sens,
           intercept,
           rwi_pred_end,
           rwi_pred_start,
           rwi_pclim_end,
           rwi_pclim_start) %>% 
    as_tibble()
  
  return(sp_fut_clim)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Compute sensitivity, RWI for each species by MC combination ---------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pull_layer <- function(brick, layer_name){
  pulled_layer <- brick %>% subset(layer_name)
}


calc_rwi_quantiles <- function(spp_code, mc_data, parallel = TRUE){
  set.seed(my_seed) # Re-setting seed at start of each iteration to ensure interrupted jobs still produce replicable results
  tic()
  
  print(paste0("Starting processing for species: ", species))
  
  ## Iterating through each species
  #print(spp_code)
  
  print("Calculating sensitivity rasters")
  
  ## Calculate n_mc versions of species' sensitivity raster
  if (parallel == TRUE) {
    sp_predictions <- mc_data %>% 
      mutate(sensitivity = future_pmap(list(sppp_code = spp_code,
                                            int_int = int_int,
                                            int_cwd = int_cwd,
                                            int_cwd2 = int_cwd2,
                                            int_pet = int_pet,
                                            int_pet2 = int_pet2,
                                            cwd_int = cwd_int,
                                            cwd_cwd = cwd_cwd,
                                            cwd_cwd2 = cwd_cwd2,
                                            cwd_pet = cwd_pet,
                                            cwd_pet2 = cwd_pet2,
                                            pet_int = pet_int,
                                            pet_cwd = pet_cwd,
                                            pet_cwd2 = pet_cwd2,
                                            pet_pet = pet_pet,
                                            pet_pet2 = pet_pet2),
                                       .f = predict_sens,
                                       .options = furrr_options(seed = my_seed, 
                                                                packages = c( "dplyr", "raster", "dtplyr"))))
  } else {
    sp_predictions <- mc_data %>% 
      mutate(sensitivity = pmap(list(sppp_code = spp_code,
                                     int_int = int_int,
                                     int_cwd = int_cwd,
                                     int_cwd2 = int_cwd2,
                                     int_pet = int_pet,
                                     int_pet2 = int_pet2,
                                     cwd_int = cwd_int,
                                     cwd_cwd = cwd_cwd,
                                     cwd_cwd2 = cwd_cwd2,
                                     cwd_pet = cwd_pet,
                                     cwd_pet2 = cwd_pet2,
                                     pet_int = pet_int,
                                     pet_cwd = pet_cwd,
                                     pet_cwd2 = pet_cwd2,
                                     pet_pet = pet_pet,
                                     pet_pet2 = pet_pet2),
                                .f = predict_sens))
  }
  
  
  print("Sensitivity rasters calculated.")
  
  sp_predictions <- sp_predictions %>% 
    mutate(cwd_const_sens = cwd_int,
           pet_const_sens = pet_int,
           int_const_sens = int_int)
  
  
  sp_predictions <- sp_predictions %>% 
    select(iter_idx, cmip_idx, sensitivity, cwd_const_sens, pet_const_sens, int_const_sens)
  gc(verbose = TRUE)
  
  print("Predicting future RWI")
  
  ## Predict future RWI for each of n_mc run
  if (parallel == TRUE){
    sp_predictions <- sp_predictions %>% 
      mutate(rwi_predictions = future_pmap(list(sppp_code = spp_code,
                                                cmip_id = cmip_idx,
                                                sensitivity = sensitivity,
                                                cwd_const_sens = cwd_const_sens,
                                                pet_const_sens = pet_const_sens,
                                                int_const_sens = int_const_sens),
                                           .f = calc_rwi_partials,
                                           .options = furrr_options(seed = my_seed,
                                                                    packages = c("raster", "dplyr", "dtplyr"))))     
  } else {
    sp_predictions <- sp_predictions %>% 
      mutate(rwi_predictions = pmap(list(sppp_code = spp_code,
                                         cmip_id = cmip_idx,
                                         sensitivity = sensitivity,
                                         cwd_const_sens = cwd_const_sens,
                                         pet_const_sens = pet_const_sens,
                                         int_const_sens = int_const_sens),
                                    .f = calc_rwi_partials)) 
  }
  
  print("Future RWI predicted.")
  
  sp_predictions <- sp_predictions %>% 
    select(iter_idx, rwi_predictions) %>% 
    unnest(rwi_predictions) %>% 
    mutate(rwi_pred_change = rwi_pred_end - rwi_pred_start,
           rwi_pclim_change = rwi_pclim_end - rwi_pclim_start, 
           rwi_pred_pclim_change_dif = rwi_pred_change - rwi_pclim_change) %>% 
    drop_na()
  
  
  gc(verbose = TRUE)
  
  print("Calculating aggregate stats")
  
  ## Calculate aggregate stats by run
  agg_stats <- sp_predictions %>% 
    # mutate(change_dif = rwi_pred_change - rwi_pclim_change) %>% 
    group_by(iter_idx) %>% 
    summarise(rwi_pred_change = mean(rwi_pred_change),
              rwi_pclim_change = mean(rwi_pclim_change))
  # change_dif = mean(change_dif))
  
  print(paste("Number of observations after summarise():", nrow(agg_stats)))
  print(paste("Number of variables after summarise():", ncol(agg_stats)))
  
  if (anyNA(agg_stats)) {
    warning("Missing or NA values found in agg_stats after calculating aggregate statistics.")
  }
  
  print("Aggregate statistics calculated.")
  
  print("Preparing historic climate data...")
  
  ## Prep historic climate data
  sp_hist <- (sp_clim %>% 
                filter(sp_code == spp_code) %>% 
                pull(clim_historic_sp))[[1]] %>% 
    rename(cwd_hist = cwd,
           pet_hist = pet)
  
  print("Historic climate data prepared")
  
  
  ## For each species, calculate cell-wise quantiles of variables from n_mc runs
  sp_predictions <- sp_predictions %>% 
    lazy_dt()
  
  print("Calculating cell-wise quantiles...")
  
  sp_predictions <- sp_predictions %>% 
    group_by(x, y) %>% 
    summarise(rwi_pred_mean = mean(rwi_pred_end),
              rwi_pred_025 = quantile(rwi_pred_end, 0.025),
              rwi_pred_050 = quantile(rwi_pred_end, 0.05),
              rwi_pred_975 = quantile(rwi_pred_end, 0.975),
              rwi_pred_950 = quantile(rwi_pred_end, 0.95),
              rwi_pclim_mean = mean(rwi_pclim_end),
              rwi_pclim_025 = quantile(rwi_pclim_end, 0.025),
              rwi_pclim_975 = quantile(rwi_pclim_end, 0.975),
              rwi_pred_change_mean = mean(rwi_pred_change),
              rwi_pred_change_025 = quantile(rwi_pred_change, 0.025),
              rwi_pred_change_050 = quantile(rwi_pred_change, 0.050),
              rwi_pred_change_950 = quantile(rwi_pred_change, 0.95),
              rwi_pred_change_975 = quantile(rwi_pred_change, 0.975),
              rwi_pclim_change_mean = mean(rwi_pclim_change),
              rwi_pclim_change_025 = quantile(rwi_pclim_change, 0.025),
              rwi_pclim_change_975 = quantile(rwi_pclim_change, 0.975),
              rwi_pred_pclim_change_dif_mean = mean(rwi_pred_pclim_change_dif),
              rwi_pred_pclim_change_dif_025 = quantile(rwi_pred_pclim_change_dif, 0.025),
              rwi_pred_pclim_change_dif_975 = quantile(rwi_pred_pclim_change_dif, 0.975),
              cwd_sens = mean(cwd_sens),
              pet_sens = mean(pet_sens),
              int_sens = mean(intercept),
              sp_code = spp_code,
              .groups = "drop")
  
  if (anyNA(sp_predictions)) {
    warning("Missing or NA values found in sp_predictions after calculating cell-wise quantiles.")
  }
  
  print("Cell-wise quantiles calculated.")
  
  ## Add back observed climate data
  sp_predictions <- sp_predictions %>% 
    left_join(sp_hist, by = c("x", "y"))
  
  if (anyNA(sp_predictions)) {
    warning("Missing or NA values found in sp_predictions after joining with historic climate data.")
  }
  
  print("Adding observed and predicted climate data")
  
  ## Add observed and predicted climate data
  sp_cmip <- (sp_clim %>% 
                filter(sp_code == spp_code) %>% 
                pull(clim_cmip_sp))[[1]]
  
  print("Added observed and predicted climate data")
  
  print("Calculating mean of cmip")
  
  sp_cmip <- sp_cmip %>%
    rowwise() %>% 
    mutate(pet_cmip_end_mean = mean(c_across(starts_with("pet_cmip_end"))),
           cwd_cmip_end_mean = mean(c_across(starts_with("cwd_cmip_end"))),
           pet_cmip_start_mean = mean(c_across(starts_with("pet_cmip_start"))),
           cwd_cmip_start_mean = mean(c_across(starts_with("cwd_cmip_start")))) %>% 
    select(x, y, cwd_cmip_start_mean, cwd_cmip_end_mean, pet_cmip_start_mean, pet_cmip_end_mean) %>% 
    drop_na()
  
   if (anyNA(sp_cmip)) {
     warning("Missing or NA values found in sp_cmip after calculating cmip means.")
   }
   
   print("Mean of cmip calculated")
   
   print("Joining cmip calculations to sp_predictions")
   
   sp_predictions <- sp_predictions %>%
     left_join(sp_cmip, by = c("x", "y")) %>% 
     as_tibble()
   
   
   print("Finished joining cmip calculations to sp_predictions")
   
   ## Write out
   write_rds(sp_predictions, paste0(output_dir, "sp_rwi_", species, ".gz"), compress = "gz")
   
   toc()
   return(agg_stats)
}

  
mc_nests <- sp_mc %>%
  group_by(sp_code) %>%
  nest() %>% 
  drop_na()


mc_nests_small <- mc_nests %>% 
  # filter(!(sp_code %in% large_range_sp)) %>% 
  mutate(predictions = pmap(list(spp_code = sp_code,
                                 mc_data = data,
                                 parallel = TRUE),
                            .f = calc_rwi_quantiles)) 

}




