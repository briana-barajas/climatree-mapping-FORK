#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Briana Barajas, Fletcher McConnell, Rosemary Juarez, Vanessa Salgado
# Project: Mapping Tree Species' Drought Sensitivity Under Climate Change
# Institution: Bren School of Environmental Science & Management - UCSB
# Date: 2024-06-07
# Purpose: Map species sensitivity and growth
#
# Input files:
# - site_pet_cwd_std_augmented (species specific)
# - sp_rwi (species specific)
# - merged_ranges_dissolve.shp
# - site_summary.csv
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#===============================================================================
# 1) Pkg imports ---------
#===============================================================================
# library(tidyverse)
# library(dbplyr)
# library(RSQLite)
# library(ggplot2)
# library(rnaturalearth)
# library(rnaturalearthdata)
# library(sf)
# library(rgeos)
# library(stringr)
# library(raster)
# library(rgdal)
# library(viridis)
# library(patchwork)
# library(effects)
# library(dplR)
# library(terra)
# select <- dplyr::select
# 
# 
base_text_size = 12
theme_set(
  theme_bw(base_size = base_text_size)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent', color=NA),
          legend.background = element_rect(fill='transparent')))

pt_size = .pt

#===============================================================================
# 2) Data imports  ---------
#===============================================================================
## Define path
# data_dir <- "~/../../capstone/climatree/raw_data/"
# output_dir <- "~/../../capstone/climatree/output/"

for(species in spp_code_list){ # <-------------------------- can choose species to run through script here
  
  # define output directory
  output_dir <- "~/../../capstone/climatree/output/intermediate-output/"
  
  # 1. Site-level regressions
  flm_df <- read_csv(paste0(output_dir, "site_pet_cwd_std_augmented_", species, ".csv")) 
  
  # 2. Species range maps
  range_file <- paste0(data_dir, 'merged_ranges_dissolve.shp')
  range_sf <- st_read(range_file)
  
  # 3. Site information
  site_smry <- read_csv(paste0(data_dir, 'site_summary.csv'))
  site_smry <- site_smry %>% 
    select(collection_id, sp_id, latitude, longitude) %>% 
    mutate(species_id = tolower(sp_id)) %>% 
    select(-sp_id)
  
  # 5. Prediction rasters
  sp_predictions <- read_rds(paste0(output_dir, "sp_rwi_", species, ".gz"))
  # rwi_list <- list.files(paste0(output_dir, "sp_rwi_pipo.gz"), pattern = ".gz", full.names = TRUE)
  # sp_predictions <- do.call('rbind', lapply(rwi_list, readRDS))
  
  
  #===============================================================================
  # Prep climate / prediction data  ---------
  #===============================================================================
  # Define species
  flm_df %>% group_by(species_id) %>% tally() %>% arrange(desc(n))
  
  #spp_code <- 'pcgl'
  spp_code <- species
  
  
  trim_df <- flm_df %>% 
    filter(species_id == spp_code)
  
  spp_predictions <- sp_predictions %>% 
    filter(sp_code == spp_code)
  
  
  #===============================================================================
  # Define example sites  ---------
  #===============================================================================
  # Pull relevant ITRDB sites
  trim_df <- trim_df %>%
    drop_na() %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
  
  high_sens = "CO559"
  low_sens = "CA585"
  
  high_coords <- trim_df %>% 
    filter(collection_id == high_sens) %>% 
    pull(geometry)
  low_coords <- trim_df %>% 
    filter(collection_id == low_sens) %>% 
    pull(geometry)
  
  high_val <- trim_df %>% 
    filter(collection_id == high_sens) %>% 
    pull(estimate_cwd.an) %>% 
    round(digits = 3)
  
  low_val <- trim_df %>% 
    filter(collection_id == low_sens) %>% 
    pull(estimate_cwd.an) %>% 
    round(digits = 3)
  
  
  high_fs <- trim_df %>% filter(collection_id == high_sens)
  low_fs <- trim_df %>% filter(collection_id == low_sens)
  
  high_lab <- paste0("sensitivity = ", as.character(high_val))
  low_lab <- paste0("sensitivity = ", as.character(low_val))
  
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
  spp_predictions <- spp_predictions %>% filter(abs(cwd_hist) < 2)
  
  # only save necessary columns for mapping
  spp_predictions <- spp_predictions %>%
    select(x, y, cwd_sens, rwi_pred_change_mean) %>% 
    mutate(species_code = species)
  
  # define directory to store final files
  output_dir <- "~/../../capstone/climatree/output/final-output/"
  
  # append the current spp_predictions to the combined data frame
  combined_predictions <- rbind(combined_predictions, spp_predictions)
  
  # save final tables in new final_output directory
  write_csv(combined_predictions, paste0(output_dir, "combined_predictions.csv"))
  write_csv(spp_predictions, paste0(output_dir, "spp_predictions_", species, ".csv"))
  
  ### Map of CWD sensitivity
  cwd_sens_map <- ggplot() +
    geom_sf(data = world) +
    geom_raster(data = spp_predictions %>% drop_na(), aes(x = x, y = y, fill = cwd_sens)) +
    theme(legend.position = c(.18,.15))+
    ylab("Latitude")+
    xlab("Longitude")+
    guides(fill=guide_legend("Sens."))+
    scale_fill_viridis(option="mako", direction = -1)+
    coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE) +
    scale_x_continuous(breaks=seq(-120,100,10)) +
    theme(axis.title.x=element_blank(),
          axis.title.y = element_blank(),
          legend.key.size = unit(8, "pt"),
          legend.title=element_text(size=base_text_size - 2), 
          legend.text=element_text(size=base_text_size - 4))+
    theme() +
    ggtitle(species)
  print(cwd_sens_map)
  
}

