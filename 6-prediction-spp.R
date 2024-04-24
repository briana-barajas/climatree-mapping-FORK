


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                               6 Predictions ----
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# clean environment
rm(list=ls())

## =======================================================================
## ---------------------------- package imports --------------------------
## =======================================================================
library(fixest)
# library(raster)
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

n_cores <- availableCores() - 12
future::plan(multisession, workers = n_cores)

my_seed <- 5597

n_mc <- 10000

## =======================================================================
## -------------------------------- load data ----------------------------
## =======================================================================
# Define path
data_dir <- "~/../../capstone/climatree/raw_data/"
output_dir <- "~/../../capstone/climatree/output/1-process-raw-data/"

# Create output directories
#out_dir <- paste0(wdir,"2_output/predictions/")
# dir.create(file.path(out_dir), showWarnings = FALSE)
#dir.create(file.path(paste0(output_dir, "sp_rwi/")), showWarnings = FALSE)
# dir.create(file.path(paste0(out_dir, "sp_hot_cells/")), showWarnings = FALSE)

# 1. Second stage model
mod_df <- read_rds(paste0(output_dir, "ss_bootstrap_pcgl.rds"))
mod_df <- mod_df %>% 
  rename(iter_idx = boot_id)

# 2. Species-standardized historic and future climate
sp_clim <- read_rds(paste0(output_dir, "sp_clim_predictions.gz")) %>% 
  filter(sp_code == "pcgl")
species_list <- sp_clim %>% select(sp_code)

## =======================================================================
## --------------------- assign MC coefs and CMIP models -----------------
## =======================================================================
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

## =======================================================================
## ------------------------- define functions ----------------------------
## =======================================================================
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
  
  ## Iterating through each species
  print(spp_code)
  
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
  sp_predictions <- sp_predictions %>% 
    mutate(cwd_const_sens = cwd_int,
           pet_const_sens = pet_int,
           int_const_sens = int_int)
  
  sp_predictions <- sp_predictions %>% 
    select(iter_idx, cmip_idx, sensitivity, cwd_const_sens, pet_const_sens, int_const_sens)
  gc(verbose = TRUE)
  
  
  
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
  sp_predictions <- sp_predictions %>% 
    select(iter_idx, rwi_predictions) %>% 
    unnest(rwi_predictions) %>% 
    mutate(rwi_pred_change = rwi_pred_end - rwi_pred_start,
           rwi_pclim_change = rwi_pclim_end - rwi_pclim_start, 
           rwi_pred_pclim_change_dif = rwi_pred_change - rwi_pclim_change)
  
  gc(verbose = TRUE)
  
  
  ## Calculate aggregate stats by run
  agg_stats <- sp_predictions %>% 
    # mutate(change_dif = rwi_pred_change - rwi_pclim_change) %>% 
    group_by(iter_idx) %>% 
    summarise(rwi_pred_change = mean(rwi_pred_change),
              rwi_pclim_change = mean(rwi_pclim_change))
  # change_dif = mean(change_dif))
  
  ## Prep historic climate data
  sp_hist <- (sp_clim %>% 
                filter(sp_code == spp_code) %>% 
                pull(clim_historic_sp))[[1]] %>% 
    rename(cwd_hist = cwd,
           pet_hist = pet)
  
  
  ## For each species, calculate cell-wise quantiles of variables from n_mc runs
  sp_predictions <- sp_predictions %>% 
    lazy_dt()
  
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
              # cwd_cmip_start = mean(cwd_cmip_start),
              # pet_cmip_start = mean(pet_cmip_start),
              # cwd_cmip_end = mean(cwd_cmip_end),
              # pet_cmip_end = mean(pet_cmip_end),
              sp_code = spp_code,
              .groups = "drop")
  
  ## Add back observed climate data
  sp_predictions <- sp_predictions %>% 
    left_join(sp_hist, by = c("x", "y"))
  
  ## Add observed and predicted climate data
  sp_cmip <- (sp_clim %>% 
                filter(sp_code == spp_code) %>% 
                pull(clim_cmip_sp))[[1]]
  
  sp_cmip <- sp_cmip %>%
    rowwise() %>% 
    mutate(pet_cmip_end_mean = mean(c_across(starts_with("pet_cmip_end"))),
           cwd_cmip_end_mean = mean(c_across(starts_with("cwd_cmip_end"))),
           pet_cmip_start_mean = mean(c_across(starts_with("pet_cmip_start"))),
           cwd_cmip_start_mean = mean(c_across(starts_with("cwd_cmip_start")))) %>% 
    select(x, y, cwd_cmip_start_mean, cwd_cmip_end_mean, pet_cmip_start_mean, pet_cmip_end_mean)
  
  sp_predictions <- sp_predictions %>% 
    left_join(sp_cmip, by = c("x", "y")) %>% 
    as_tibble()
  
  ## Write out
  sp_predictions %>% 
    write_rds(paste0(output_dir, "sp_rwi_pcgl.gz"), compress = "gz")
  
  toc()
  return(agg_stats)
}