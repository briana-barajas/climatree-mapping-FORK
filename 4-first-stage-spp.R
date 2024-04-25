5A

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 5/17/20
# Purpose: Run regressions to explore impact of historical climate on weather sensitivity
#
# Input files:
# - site_pet_cwd_std.csv: Table of first stage regression parameters for baseline specification. Generated in "4a. First stage.R"
# - site_ave_clim.gz: Site-level climate parameters. Generated using "3b. Species_niche.R"
# - site_summary.csv: Summary data about each site. Generated using "1b. Parse ITRDB.R"
# - species_gen_gr.csv: Lookup key for angiosperm / gymnosperm differentiation
# 
# Output files:
# - site_pet_cwd_std_augmented.csv: Intermediate file to enable downstream identification of outlier assignment.
# - mc_sample.gz: Intermediate file of all bootstrapped samples.
# - ss_bootstrap.rds: Final second stage model parameters for each of n_mc Monte Carlo runs.
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(MASS)
library(tidyverse)
library(broom)
library(purrr)
library(margins)
library(tidylog)
library(fixest)
library(gstat)
library(sf)
library(units)
library(dtplyr)
library(marginaleffects)

set.seed(5597)

select <- dplyr::select

n_mc <- 10000


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Define path
data_dir <- "~/../../capstone/climatree/raw_data/"
output_dir <- "~/../../capstone/climatree/output/1-process-raw-data/"

# 1. Site-level regressions
flm_df <- read_csv(paste0(output_dir, 'site_pet_cwd_std.csv'))

# 2. Historic site-level climate
ave_site_clim <- read_rds(paste0(output_dir, "site_ave_clim.gz"))
flm_df <- flm_df %>% 
  left_join(ave_site_clim, by = c("collection_id"))

# 3. Site information
site_df <- read_csv(paste0(data_dir, 'site_summary.csv'))
site_df <- site_df %>% 
  select(collection_id, sp_id, latitude, longitude)
site_df <- site_df %>% 
  rename(species_id = sp_id) %>% 
  mutate(species_id = str_to_lower(species_id))

# 4. Species information
sp_info <- read_csv(paste0(data_dir, 'species_metadata.csv'))
sp_info <- sp_info %>% 
  select(species_id, genus, gymno_angio, family)
site_df <- site_df %>% 
  left_join(sp_info, by = "species_id")

# Merge back into main flm_df
flm_df <- flm_df %>% 
  left_join(site_df, by = "collection_id") %>% 
  filter(species_id == "pcgl") #<<<<<<<<<<<<<<<<<<<------------ input desired species code(s) here to run script------------------------------

# define new folder to store outputs
output_dir <- "~/../../capstone/climatree/output/test-output/"

# re-define flm_df for for loop
original_flm_df <- flm_df

# define species_id column to iterate through for for loop
spp_code_list <- flm_df$species_id

# for loop for iterating script through all species

for(i in spp_code_list) {
  
  
  flm_df <- original_flm_df %>% 
    filter(species_id == i) # filter for only rows with species of interest
  
  if (nrow(flm_df) == 0) {
    print(paste0("No data found for species ", i))
    next
  }
  
  #Filter for species code pcgl (White Spruce)
  #flm_df <- flm_df %>% 
  #filter(species_id == "psme")
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Prep and trim data -----------------------------------------------------
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Add weighting based on inverse of first stage variance
  flm_df <- flm_df %>% 
    mutate(cwd_errorweights = 1 / (std.error_cwd.an),
           errorweights2 = sqrt(ntrees),
           pet_errorweights = 1 / (std.error_pet.an),
           int_errorweights = 1 / (std.error_intercept))
  
  # Identify and trim extreme outliers
  cwd_est_bounds = quantile(flm_df$estimate_cwd.an, c(0.01, 0.99),na.rm=T)
  pet_est_bounds = quantile(flm_df$estimate_pet.an, c(0.01, 0.99),na.rm=T)
  cwd_spstd_bounds = quantile(flm_df$cwd.spstd, c(0.01, 0.99), na.rm = T)
  pet_spstd_bounds = quantile(flm_df$pet.spstd, c(0.01, 0.99), na.rm = T)
  
  
  
  flm_df <- flm_df %>%
    mutate(outlier = (estimate_cwd.an<cwd_est_bounds[1]) |
             (estimate_cwd.an>cwd_est_bounds[2]) |
             (estimate_pet.an<pet_est_bounds[1]) |
             (estimate_pet.an>pet_est_bounds[2]))
  
  # Save out full flm_df to simplify downstream scripts and ensure consistency
  flm_df %>% write.csv(paste0(output_dir, "site_pet_cwd_std_augmented_", i, ".csv"))
  #flm_df %>% write.csv(paste0(output_dir, "site_pet_cwd_std_augmented_pcgl.csv"))
  
  # Trim outliers
  trim_df <- flm_df %>% 
    filter(outlier==0) %>% 
    drop_na()
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Spatial autocorrelation of trim_df ---------------------------------
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  site_points=st_as_sf(trim_df,coords=c("longitude","latitude"),crs="+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
  
  vg <-variogram(estimate_cwd.an~1, site_points, cutoff = 1500, width = 10)
  vg.fit <- fit.variogram(vg, model = vgm(1, "Sph", 900, 1))
  plot(vg, vg.fit)
  # print(paste0("Range before hitting sill (km): "), as.character(vg.fit[2,3]))
  
  vg.range = vg.fit[2,3] * 1000
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Quick test of primary regression ---------------------------------
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  formula = as.formula("estimate_cwd.an ~ cwd.spstd + (cwd.spstd^2) + pet.spstd + (pet.spstd^2)")
  mod_data <- trim_df
  cwd_mod <- feols(formula, data = mod_data, weights = mod_data$cwd_errorweights,
                   vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
  summary(cwd_mod)
  
  marg_fx_df <- function(mod){
    inc <- 0.1
    min <- -2.5
    max <- 2.5
    cwd_pred <- predictions(mod, newdata = datagrid(pet.spstd = 0, cwd.spstd = seq(min,max,inc))) %>% 
      mutate(variation = "cwd")
    pet_pred <- predictions(mod, newdata = datagrid(pet.spstd = seq(min,max,inc), cwd.spstd = 0)) %>% 
      mutate(variation = "pet")
    return(rbind(cwd_pred, pet_pred))
  }
  
  
  preds <- marg_fx_df(cwd_mod)
  
  cwd_mfx_plot <- preds %>% 
    filter(variation == "cwd") %>% 
    ggplot(aes(x = cwd.spstd)) + 
    geom_line(aes(y = estimate)) +
    geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2)
  cwd_mfx_plot
  
  pet_mfx_plot <- preds %>% 
    filter(variation == "pet") %>% 
    ggplot(aes(x = pet.spstd)) + 
    geom_line(aes(y = estimate)) +
    geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2)
  pet_mfx_plot
  
  formula = as.formula("estimate_pet.an ~ cwd.spstd + pet.spstd + (cwd.spstd^2) + (pet.spstd^2)")
  pet_mod <- feols(formula, weights = mod_data$pet_errorweights, data = mod_data,
                   vcov = conley(cutoff = vg.range/1000, distance = "spherical"))
  summary(pet_mod)
  preds <- marg_fx_df(pet_mod)
  
  cwd_mfx_plot <- preds %>% 
    filter(variation == "cwd") %>% 
    ggplot(aes(x = cwd.spstd)) + 
    geom_line(aes(y = estimate)) +
    geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2)
  cwd_mfx_plot
  
  pet_mfx_plot <- preds %>% 
    filter(variation == "pet") %>% 
    ggplot(aes(x = pet.spstd)) + 
    geom_line(aes(y = estimate)) +
    geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2)
  pet_mfx_plot
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Identify spatially proximate blocks of sites ---------------------------------
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  site_list <- trim_df %>%
    pull(collection_id) %>%
    unique()
  n_sites <- length(site_list)
  
  site_dist=st_distance(site_points)
  rownames(site_dist)=site_points$collection_id
  colnames(site_dist)=site_points$collection_id
  # save(site_dist,file=paste0(wdir,"out/site_distances.Rdat"))
  # load(paste0(wdir,"out/site_distances.Rdat"))
  
  dist_df <- as_tibble(site_dist) %>% 
    drop_units() 
  
  dist_df <- dist_df %>%
    lazy_dt() %>% 
    mutate(collection_id = names(dist_df)) %>% 
    # select(collection_id, site_list) %>% 
    # filter(collection_id %in% site_list) %>% 
    mutate(across(.cols = !collection_id, ~(.x < vg.range))) %>% 
    # mutate(across(.cols = !collection_id, ~ifelse((.x < range), collection_id, "DROP"))) %>% 
    as_tibble()
  
  block_list <- c()
  for (site in site_list){
    block_sites <- dist_df %>% 
      filter(get(site) == TRUE) %>% 
      pull(collection_id)
    block_list[site] <- list(block_sites)
  }
  #save(block_list,file=paste0(output_dir,"spatial_blocks_pcgl.Rdat"))
  # load(file=paste0(wdir,"out/second_stage/spatial_blocks.Rdat"))
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Create block bootstrap draws  ---------------------------------
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  draw_blocks <- function(site_sample){
    samp <- site_sample %>% pull(samp)
    blocked_draw <- (block_list[samp] %>% unlist() %>% unname())[1:n_sites]
    return(blocked_draw)  
  }
  
  n_obs = n_mc * n_sites
  block_draw_df <- tibble(boot_id = rep(1:n_mc, each = n_sites)) %>% 
    lazy_dt() %>% 
    mutate(samp = sample(site_list,size=n_obs,replace=TRUE)) %>%
    group_by(boot_id) %>% 
    nest()
  
  block_draw_df <- block_draw_df %>% 
    mutate(sites = map(.x = data, .f = draw_blocks)) %>% # COULD PARALLELIZE HERE?
    select(boot_id, sites) %>% 
    as_tibble() %>% 
    unnest(sites) 
  
  block_draw_df <- block_draw_df %>% 
    rename(collection_id = sites)
  
  ## Identify number of draws needed for each site
  n_draws <- block_draw_df %>% 
    group_by(collection_id) %>% 
    tally() %>% 
    rename(n_draw = n)
  
  trim_df <- trim_df %>% 
    left_join(n_draws, by = "collection_id")
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Random draws of coefs from first stage ---------------------------------
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## Function to create random draws of first stage coefficients
  draw_coefs <- function(n, cwd_est, pet_est, int_est, cwd_ste, pet_ste, int_ste, 
                         cwd_pet_cov, pet_int_cov, cwd_int_cov){
    mu <- c("cwd_coef" = cwd_est, 
            "pet_coef" = pet_est, 
            "int_coef" = int_est)
    vcov <- matrix(data = NA, nrow = 3, ncol = 3)
    vcov[1,1] <- cwd_ste^2
    vcov[2,2] <- pet_ste^2
    vcov[3,3] <- int_ste^2
    vcov[1,2] <- cwd_pet_cov
    vcov[2,1] <- cwd_pet_cov
    vcov[1,3] <- cwd_int_cov
    vcov[3,1] <- cwd_int_cov
    vcov[2,3] <- pet_int_cov
    vcov[3,2] <- pet_int_cov
    
    draw <- mvrnorm(n, mu, vcov)
    draw <- as_tibble(draw)
    draw$iter_idx <- seq(1,n)
    draw <- draw %>% select(iter_idx, cwd_coef, pet_coef, int_coef)
    return(draw)
  }
  
  ## Create needed number (n_draw) of random draws of first stage coefficients for each site
  trim_df <- trim_df %>% 
    drop_na()
  
  mc_df <- trim_df %>%
    mutate(coef_draws = pmap(list(n = trim_df$n_draw + 1, 
                                  cwd_est = trim_df$estimate_cwd.an, 
                                  pet_est = trim_df$estimate_pet.an,
                                  int_est = trim_df$estimate_intercept, 
                                  cwd_ste = trim_df$std.error_cwd.an,
                                  pet_ste = trim_df$std.error_pet.an, 
                                  int_ste = trim_df$std.error_intercept,
                                  cwd_pet_cov = trim_df$cov_cwd_pet, 
                                  cwd_int_cov = trim_df$cov_int_cwd,
                                  pet_int_cov = trim_df$cov_int_pet), 
                             draw_coefs))
  
  
  ## Unnest to create dataframe of n_site X n_draw coefficient estimates
  mc_df <- mc_df %>% 
    unnest(coef_draws) %>% 
    select(collection_id, iter_idx, cwd_coef, pet_coef, int_coef, cwd.spstd, 
           pet.spstd, latitude, longitude, cwd_errorweights, pet_errorweights, int_errorweights)
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Merge first stage draws back to bootstrap dataframe -------------------------
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  block_draw_df <- block_draw_df %>% 
    group_by(collection_id) %>% 
    mutate(iter_idx = 1:n())
  
  block_draw_df <- block_draw_df %>% 
    left_join(mc_df, by = c("collection_id", "iter_idx"))
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Export first stage draws to pull summary stats -------------------------
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #block_draw_df %>% 
  # select(boot_id, collection_id, cwd_coef, pet_coef, int_coef, cwd.spstd, pet.spstd) %>% 
  #write_rds(paste0(output_dir, "mc_sample_pcgl.gz"), compress = "gz")
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Run bootstrap estimation of second stage model -------------------------
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  
  bs_constant <- function(data){
    # data <- data %>% # Needed to add this since block bootstrap is returning nested tibble
    #   unnest(cols = c(data))
    const_sens <- data %>% 
      summarise(cwd_const_sens = weighted.mean(cwd_coef, cwd_errorweights),
                pet_const_sens = weighted.mean(pet_coef, pet_errorweights),
                int_const_sens = weighted.mean(int_coef, pet_errorweights))
    return(const_sens)
  }
  
  
  ## Create dataframe holding bootstrap samples
  boot_df <- block_draw_df %>%
    select(-iter_idx) %>% 
    group_by(boot_id) %>% 
    nest()
  
  
  ## Estimate second stage models
  boot_df <- boot_df %>% 
    mutate(const_sens = map(.x = data, .f = bs_constant)) %>%
    unnest_wider(const_sens) %>%
    mutate(estimates = map(.x = data, .f = bs_ss)) %>% 
    unnest_wider(estimates) %>% 
    select(-data) %>% 
    ungroup()
  
  
  ## Save out bootstrapped coefficients
  write_rds(boot_df, paste0(output_dir, "ss_bootstrap_", i, ".rds"))
  #write_rds(boot_df, paste0(output_dir, "ss_bootstrap_pcgl.rds"))
  
  print(paste0("Processing of tree species ", i, " is complete."))
  
}






6 
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
6 





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Joan Dudney, Frances Moore
# Project: Treeconomics
# Date: 5/27/20
# Purpose: Create predictions of growth impacts from climate change
#
# Input files:
# - ss_bootstrap.rds: R model object saved from Second stage
# - sp_clim_predictions.gz:
# 
# Output files
# - sp_rwi/<sp_code>.gz:
#
# ToDo:
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
data_dir <- "~/../../capstone/climatree/raw_data/"
output_dir1 <- "~/../../capstone/climatree/output/1-process-raw-data/"
output_dir2 <- "~/../../capstone/climatree/output/test-output/"

# Create output directories
#out_dir <- paste0(wdir,"2_output/predictions/")
# dir.create(file.path(out_dir), showWarnings = FALSE)
#dir.create(file.path(paste0(output_dir, "sp_rwi/")), showWarnings = FALSE)
# dir.create(file.path(paste0(out_dir, "sp_hot_cells/")), showWarnings = FALSE)

for (i in spp_code_list) {
  
  
  # 1. Second stage model
  mod_df <- read_rds(paste0(output_dir2, "ss_bootstrap_", i, ".rds"))
  mod_df <- mod_df %>% 
    rename(iter_idx = boot_id)
  
  # 2. Species-standardized historic and future climate
  sp_clim <- read_rds(paste0(output_dir1, "sp_clim_predictions.gz")) %>% 
    filter(sp_code == i)
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
      write_rds(paste0(output_dir2, "sp_rwi_", i, ".gz"), compress = "gz")
    
   # toc()
    #return(agg_stats)
  }
  
}

#mc_nests <- sp_mc %>%
# group_by(sp_code) %>%
#nest() %>% 
#drop_na()

#mc_nests_small <- mc_nests %>% 
# filter(!(sp_code %in% large_range_sp)) %>% 
# mutate(predictions = pmap(list(spp_code = sp_code,
#                               mc_data = data,
#                              parallel = TRUE),
#                        .f = calc_rwi_quantiles)) 

#agg_stats <- mc_nests_small %>% 
# select(-data) %>% 
#unnest(predictions) %>% 
#write_rds(file = paste0(output_dir, "mc_agg_stats_pcgl.gz"), compress = "gz")


#test <- agg_stats %>% 
# group_by(iter_idx) %>% 
#summarise(rwi_pred_change = mean(rwi_pred_change))
#test %>%
# pull(rwi_pred_change) %>% 
#quantile(c(0.025, 0.5, 0.975))




7
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
7





#===============================================================================
# Species deep dive plots
# 
# 
#===============================================================================

#===============================================================================
# 1) Pkg imports ---------
#===============================================================================
library(tidyverse)
library(dbplyr)
library(RSQLite)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(rgeos)
library(stringr)
library(raster)
library(rgdal)
library(viridis)
library(patchwork)
library(effects)
library(dplR)
library(terra)
select <- dplyr::select


base_text_size = 12
theme_set(
  theme_bw(base_size = base_text_size)+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          # text=element_text(family ="Helvetica"),
          panel.background = element_rect(fill='transparent'), 
          plot.background = element_rect(fill='transparent', color=NA), 
          legend.background = element_rect(fill='transparent')))

pt_size = .pt

#===============================================================================
# 2) Data imports  ---------
#===============================================================================
### Define path
data_dir <- "~/../../capstone/climatree/raw_data/"
output_dir <- "~/../../capstone/climatree/output/1-process-raw-data/"


for(i in spp_code_list) {
  
  # 1. Site-level regressions
  flm_df <- read_csv(paste0(output_dir, "site_pet_cwd_std_augmented_", i, ".csv")) 
  
  # 2. Species range maps
  range_file <- paste0(data_dir, 'merged_ranges_dissolve.shp')
  range_sf <- st_read(range_file)
  
  # 3. Site information
  site_smry <- read_csv(paste0(data_dir, 'site_summary.csv'))
  species_metadata <- read_csv(paste0(data_dir, 'species_metadata.csv'))
  
  site_smry <- site_smry %>% 
    select(collection_id, sp_id, latitude, longitude) %>% 
    mutate(species_id = tolower(sp_id)) %>% 
    select(-sp_id)
  # site_loc <- site_smry %>% 
  #   select(collection_id, latitude, longitude)
  # flm_df <- flm_df %>% 
  #   left_join(site_loc, by = "collection_id")
  
  # # 4. Species information
  # sp_info <- read_csv(paste0(wdir, 'species_gen_gr.csv'))
  # sp_info <- sp_info %>% 
  #   select(species_id, genus, gymno_angio, family)
  # site_smry <- site_smry %>% 
  #   left_join(sp_info, by = c("species_id"))
  
  # 5. Prediction rasters
  #rwi_list <- list.files(paste0(output_dir, "sp_rwi/"), pattern = ".gz", full.names = TRUE)
  #sp_predictions <- do.call('rbind', lapply(rwi_list, readRDS))
  sp_predictions <- read_rds(paste0(output_dir, "sp_rwi_", i, ".gz"))
  
  # 6. Dendro examples - note: exporting two pipo sites in first stage script
  dendro_ex <- read_csv(paste0(output_dir, "example_sites.csv"))
  
  # 7. Raw dendro file for one site
  #rwl_path <- paste0(data_dir, "ca585.rwl")
  #rwl_dat <- read.tucson(paste0(rwl_path))
  
  
  
  #===============================================================================
  # Prep climate / prediction data  ---------
  #===============================================================================
  # Define species
  flm_df %>% group_by(species_id) %>% tally() %>% arrange(desc(n))
  
  #spp_code <- 'pcgl'
  spp_code <- i
  
  
  trim_df <- flm_df %>% 
    filter(species_id == spp_code)
  
  spp_predictions <- sp_predictions %>% 
    filter(sp_code == spp_code)
  
  # sp_fut <- (spp_predictions %>% 
  #              pull(clim_future_sp))[[1]]
  # names(sp_fut) <- c("cwd.fut", "pet.fut")
  # 
  # sp_hist <- (spp_predictions %>% 
  #               pull(clim_historic_sp))[[1]]
  # sp_sens  <- (spp_predictions %>% 
  #                pull(sensitivity))[[1]]
  # sp_rwi  <- (spp_predictions %>% 
  #               pull(rwi_predictions))[[1]]
  # names(sp_rwi) <- "rwi"
  # 
  # 
  # clim_compare <- brick(c(sp_fut, sp_hist, sp_sens, sp_rwi))
  # clim_compare <- clim_compare %>% 
  #   as.data.frame(xy = TRUE)
  
  
  #===============================================================================
  # Define example sites  ---------
  #===============================================================================
  # Pull relevant ITRDB sites
  #trim_df <- trim_df %>%
  #drop_na() %>%
  #st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
  
  #high_sens = "CO559"
  #low_sens = "CA585"
  
  #high_coords <- trim_df %>% 
  #filter(collection_id == high_sens) %>% 
  #pull(geometry)
  #low_coords <- trim_df %>% 
  #filter(collection_id == low_sens) %>% 
  #pull(geometry)
  
  #high_val <- trim_df %>% 
  #filter(collection_id == high_sens) %>% 
  #pull(estimate_cwd.an) %>% 
  #round(digits = 3)
  
  #low_val <- trim_df %>% 
  #filter(collection_id == low_sens) %>% 
  #pull(estimate_cwd.an) %>% 
  #round(digits = 3)
  
  
  #high_fs <- trim_df %>% filter(collection_id == high_sens)
  #low_fs <- trim_df %>% filter(collection_id == low_sens)
  
  #high_lab <- paste0("sensitivity = ", as.character(high_val))
  #low_lab <- paste0("sensitivity = ", as.character(low_val))
  
  #high_color <- "#404788"
  #low_color <- "#efca2a"
  #low_color <- "#1b9e77"
  
  #440154FF
  
  #===============================================================================
  # Step 1: data and detrending  ---------
  #===============================================================================
  ### Map of ITRDB sites and species range
  
  # Pull relevant range map
  sp_range <- range_sf %>% 
    filter(sp_code == spp_code)
  sp_bbox <- st_bbox(sp_range)
  
  lon_lims <- c(sp_bbox$xmin - 1, sp_bbox$xmax + 1)
  lat_lims <- c(sp_bbox$ymin - 1, sp_bbox$ymax + 1)
  
  # Plot species ranges
  world <- ne_coastline(scale = "medium", returnclass = "sf")
  map <- ggplot(trim_df, aes(x = Longitude, y = Latitude))
  
  #===============================================================================
  # Step 5: Prediction of sensitivity  ---------
  #===============================================================================
  spp_predictions <- spp_predictions %>% filter(abs(cwd_hist) < 2)  ## TODO - Figure out correct cut-off for predictions
  
  spp_predictions <- spp_predictions %>%
    select(x, y, cwd_sens, rwi_pred_change_mean)
  
  # create rasters for each species to save
  # v <- vect(spp_predictions, geom = c("x", "y"), crs = "EPSG:4326")
  # r <- rast(v, res = 0.5) 
  #rasterized_raster <- rasterize(v, r, field = "cwd_sens")
  
  
  # specify new directory to store final table values for each species
  output_dir2 <- "~/../../capstone/climatree/output/test-output/"
  
  # save final tables in new final_output directory
  write_rds(spp_predictions, paste0(output_dir2, "spp_predictions_", i, ".rds"))
  
  ### Map of CWD sensitivity
  cwd_sens_map <- ggplot() +
    geom_sf(data = world) +
    geom_raster(data = spp_predictions %>% drop_na(), aes(x = x, y = y, fill = cwd_sens)) +
    #theme_bw(base_size = 22)+
    theme(legend.position = c(.18,.15))+
    ylab("Latitude")+
    xlab("Longitude")+
    guides(fill=guide_legend("Sens."))+
    #scale_fill_viridis_c(direction = -1) +
    scale_fill_viridis(option="mako", direction = -1)+
    coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE) +
    scale_x_continuous(breaks=seq(-120,100,10)) +
    theme(axis.title.x=element_blank(),
          axis.title.y = element_blank(),
          # axis.text.x=element_text(size=base_text_size - 6),
          # axis.text.y=element_text(size = base_text_size - 6),
          legend.key.size = unit(8, "pt"),
          legend.title=element_text(size=base_text_size - 2), 
          legend.text=element_text(size=base_text_size - 4))+
    theme()
  cwd_sens_map
  
  #===============================================================================
  # Step 5: Prediction of RWI change  ---------
  #===============================================================================
  ### Map of CWD change
  #spp_predictions <- spp_predictions %>% 
  #mutate(cwd_change = cwd_cmip_end_mean - cwd_cmip_start_mean,
  #pet_change = pet_cmip_end_mean - pet_cmip_start_mean)
  
  
  #cwd_change_map <- ggplot() +
  #geom_sf(data = world) +
  #geom_raster(data = spp_predictions %>% drop_na(), aes(x = x, y = y, fill = cwd_change)) +
  #theme_bw(base_size = 22)+
  #guides(fill=guide_legend("Δ CWD"))+
  #theme(legend.position = c(.18,.15))+
  #ylab("Latitude")+
  #xlab("Longitude")+
  #scale_fill_viridis_c(direction = -1) +
  #scale_fill_viridis(option="magma")+
  #coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE) +
  #scale_x_continuous(breaks=seq(-120,100,10)) +
  #theme(axis.title.x=element_blank(),
  #axis.title.y = element_blank(),
  # axis.text.x=element_text(size=base_text_size - 6),
  # axis.text.y=element_text(size = base_text_size - 6),
  #legend.key.size = unit(8, "pt"),
  #legend.title=element_text(size=base_text_size - 2), 
  #legend.text=element_text(size=base_text_size - 4))
  #cwd_change_map
  
  
  ### Map of predicted RWI
  rwi_map <- ggplot() +
    geom_sf(data = world) +
    geom_raster(data = spp_predictions %>% drop_na(), aes(x = x, y = y, fill = rwi_pred_change_mean)) +
    # theme_bw(base_size = 12)+
    ylab("Latitude")+
    xlab("Longitude")+
    scale_fill_viridis_c(direction = -1) +
    #scale_fill_viridis(option="mako")+
    guides(fill=guide_legend(title="Δ RWI"))+
    coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE) +
    scale_x_continuous(breaks=seq(-120,100,10)) +
    theme(legend.position = c(.18,.15),
          # axis.text.x=element_text(size=base_text_size - 6),
          # axis.text.y=element_text(size = base_text_size - 6),
          axis.title.x=element_blank(),
          axis.title.y = element_blank(),
          legend.key.size = unit(8, "pt"),
          legend.title=element_text(size=base_text_size - 2), 
          legend.text=element_text(size=base_text_size - 4))
  rwi_map
  
}