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
          # text=element_text(family ="Helvetica"),
          panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent', color=NA),
          legend.background = element_rect(fill='transparent')))

pt_size = .pt

#===============================================================================
# 2) Data imports  ---------
#===============================================================================
### Define path
# data_dir <- "~/../../capstone/climatree/raw_data/"
# input_dir <- "~/../../capstone/climatree/output/new-output/"
#output_dir <- "~/../../capstone/climatree/output/intermediate-output/"

for(species in spp_code_list){
  
  # define output directory
  output_dir <- "~/../../capstone/climatree/output/intermediate-output/"
  
# 1. Site-level regressions
flm_df <- read_csv(paste0(output_dir, "site_pet_cwd_std_augmented_", species, ".csv")) 

# 2. Species range maps
range_file <- paste0(input_dir, 'merged_ranges_dissolve.shp')
range_sf <- st_read(range_file)

# 3. Site information
site_smry <- read_csv(paste0(input_dir, 'site_summary.csv'))
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
sp_predictions <- read_rds(paste0(output_dir, "sp_rwi_", species, ".gz"))
# rwi_list <- list.files(paste0(output_dir, "sp_rwi_pipo.gz"), pattern = ".gz", full.names = TRUE)
# sp_predictions <- do.call('rbind', lapply(rwi_list, readRDS))

# 6. Dendro examples - note: exporting two pipo sites in first stage script
#dendro_ex <- read_csv(paste0(output_dir, "example_sites.csv"))

# 7. Raw dendro file for one site
#rwl_path <- paste0(data_dir, "ca585.rwl")
#rwl_dat <- read.tucson(paste0(rwl_path))



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

high_color <- "#404788"
  low_color <- "#efca2a"
    low_color <- "#1b9e77"
      
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
  theme() +
  ggtitle(species)
print(cwd_sens_map)

}
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
# rwi_map <- ggplot() +
#   geom_sf(data = world) +
#   geom_raster(data = spp_predictions %>% drop_na(), aes(x = x, y = y, fill = rwi_pred_change_mean)) +
#   # theme_bw(base_size = 12)+
#   ylab("Latitude")+
#   xlab("Longitude")+
#   scale_fill_viridis_c(direction = -1) +
#   #scale_fill_viridis(option="mako")+
#   guides(fill=guide_legend(title="Δ RWI"))+
#   coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE) +
#   scale_x_continuous(breaks=seq(-120,100,10)) +
#   theme(legend.position = c(.18,.15),
#         # axis.text.x=element_text(size=base_text_size - 6),
#         # axis.text.y=element_text(size = base_text_size - 6),
#         axis.title.x=element_blank(),
#         axis.title.y = element_blank(),
#         legend.key.size = unit(8, "pt"),
#         legend.title=element_text(size=base_text_size - 2), 
#         legend.text=element_text(size=base_text_size - 4)) +
#   ggtitle(species)
# print(rwi_map)
    
