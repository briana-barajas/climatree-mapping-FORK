library(tidyverse)
library(here)

library(sf)     # vector data
library(terra)  # raster data
library(tmap)   # mapping
#--------------------------------------------------------------------------
#to map out the maps for our projects and presentation.

#more so on the top 26 species

meta_df <- read_csv("~/../../capstone/climatree/input/external/species_metadata.csv") %>% 
  rename(sp_code = species_id) %>% 
  dplyr::select(sp_code, spp, common_name)

combined_pred_26 <- read_csv("~/../../capstone/climatree/output/final-output/combined_predictions.csv") %>% 
  rename(sp_code = species_code) %>% 
  left_join(meta_df, join_by(sp_code)) %>% 
  select(-x, -y, -cwd_sens) %>% 
  distinct(sp_code, .keep_all = TRUE)

#----------------------------------------------------------------------------------
all_species_range <- st_read(here("~", "..", "..", "capstone", "climatree", "input", "external",  "merged_ranges_dissolve.shp")) %>% 
  st_make_valid()


top_26_species_range <- all_species_range %>%
  filter(sp_code %in% combined_pred_26$sp_code)

top_range <- st_union(top_26_species_range)

merged_sf <- st_sf(geometry = st_sfc(top_range))

data(World)

range_26 <- tm_shape(World) +
  tm_fill("#ffffdd") +
  tm_borders() +
  # tm_polygons() +
  # tm_style("classic") +
  tm_shape(top_26_species_range)+
  tm_fill(col = '#364030', border.alpha = 0 )

range_26

tmap_save(range_26, 'range_26.jpg')

