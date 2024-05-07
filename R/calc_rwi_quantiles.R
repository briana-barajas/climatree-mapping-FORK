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
  
  print(paste("Number of observations after summarise():", nrow(agg_stats)))
  print(paste("Number of variables after summarise():", ncol(agg_stats)))
  
  ## Prep historic climate data
  sp_hist <- (sp_clim %>% 
                filter(sp_code == spp_code) %>% 
                pull(clim_historic_sp))[[1]] %>% 
    rename(cwd_hist = cwd,
           pet_hist = pet)
  
  assign("sp_hist", sp_hist, envir = .GlobalEnv)
  
  # ## Write out full mc rwi change results for subset of hot cells (pet ~= 1)
  # hot_cells <- sp_hist %>% filter(pet_hist > 0.9, pet_hist < 1.1)
  # hot_cells <- hot_cells %>% 
  #   lazy_dt() %>% 
  #   left_join(sp_predictions, by = c("x", "y")) %>% 
  #   mutate(sp_code = spp_code) %>% 
  #   select(sp_code, iter_idx, x, y, cwd_hist, pet_hist, rwi_pred_change) %>% 
  #   as.data.frame()
  # hot_cells %>% 
  #   write_rds(file = paste0(out_dir, "sp_hot_cells/", spp_code, ".gz"), compress = "gz")
  # 
  # ## Contrast RWI change in wettest and dryest sites (all warm)
  # pet_range <- sp_hist %>% pull(pet_hist) %>% quantile(c(0.75, 1))
  # hot_cells <- sp_hist %>% filter(pet_hist > pet_range[1], pet_hist < pet_range[2])
  # cwd_quantile <- hot_cells %>% pull(cwd_hist) %>% quantile(c(0.1, 0.9))
  # wet_cells <- hot_cells %>% filter(cwd_hist < cwd_quantile[1]) %>% 
  #   mutate(wet_dry = "wet")
  # dry_cells <- hot_cells %>% filter(cwd_hist > cwd_quantile[2]) %>% 
  #   mutate(wet_dry = "dry")
  # extreme_cells <- rbind(wet_cells, dry_cells)
  # 
  # extreme_cells <- extreme_cells %>% 
  #   lazy_dt() %>% 
  #   left_join(sp_predictions, by = c("x", "y")) %>% 
  #   group_by(iter_idx, wet_dry) %>% 
  #   summarise(rwi_change = mean(rwi_pred_change)) %>% 
  #   pivot_wider(names_from = wet_dry, values_from = rwi_change) %>% 
  #   mutate(wet_dry_dif = wet - dry) %>% 
  #   as.data.frame()
  #   
  # agg_stats <- agg_stats %>%
  #   left_join(extreme_cells %>% select(iter_idx, wet_dry_dif), by = "iter_idx")
  
  
  # ## Contrast RWI change under two scenarios for PET-centered sites
  # pet_range = sp_hist$pet_hist %>% quantile(c(0.45, 0.55))
  # pet_cells <- sp_hist %>% filter(pet_hist > pet_range[1], pet_hist < pet_range[2])
  # cwd_quantile <- pet_cells %>% pull(cwd_hist) %>% quantile(c(0.1, 0.9))
  # wet_cells <- pet_cells %>% filter(cwd_hist < cwd_quantile[1]) %>% 
  #   mutate(wet_dry = "wet")
  # dry_cells <- pet_cells %>% filter(cwd_hist > cwd_quantile[2]) %>% 
  #   mutate(wet_dry = "dry")
  # extreme_cells <- rbind(wet_cells, dry_cells)
  # 
  # extreme_cells <- extreme_cells %>% 
  #   lazy_dt() %>% 
  #   left_join(sp_predictions, by = c("x", "y")) %>% 
  #   group_by(iter_idx, wet_dry) %>% 
  #   summarise(rwi_pred_change = mean(rwi_pred_change),
  #             rwi_pclim_change = mean(rwi_pclim_change)) %>% 
  #   pivot_wider(names_from = wet_dry, values_from = c(rwi_pclim_change, rwi_pred_change)) %>% 
  #   mutate(wet_pred_pclim_dif = rwi_pred_change_wet - rwi_pclim_change_wet,
  #          dry_pred_pclim_dif = rwi_pred_change_dry - rwi_pclim_change_dry) %>% 
  #   as.data.frame()
  # 
  # agg_stats <- agg_stats %>%
  #   left_join(extreme_cells %>% select(iter_idx, wet_pred_pclim_dif, dry_pred_pclim_dif),
  #             by = "iter_idx")
  
  
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
  
  assign("sp_cmip", sp_cmip, envir = .GlobalEnv)
  
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
  write_rds(sp_predictions, paste0(output_dir, "sp_rwi_pipo_old.gz"), compress = "gz")
  
  toc()
  return(agg_stats)
}