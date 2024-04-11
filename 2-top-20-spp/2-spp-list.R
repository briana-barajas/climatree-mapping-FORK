#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Rosemary Juarez, Briana Barajas, Fletcher McConnell, Vanessa Salgado
# Project: Mapping Global Tree Vulnerability Under Climate Change
# Date: 2024-01-28
# Purpose: 1) Isolate the 20 most common tree species in the ITRDB
#
# Input files:
#   site_summary.csv
#   rwi_long.csv
#   species_metadata.csv
# 
# Output files:
#   NA
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## =======================================================================
## -------------------------- Package imports ----------------------------
## =======================================================================

library(tidyverse)
library(here)
library(gghighlight)

## =======================================================================
## ---------------------------- Load data --------------------------------
## =======================================================================

# set working directory
remote_wdir <- here("~","..", "..", "capstone", "climatree", "raw_data")

# read in tree ring data
rwi_long <- read_csv(here(remote_wdir, "rwi_long.csv"))

# read in site summary table
site_smry <- read_csv(here(remote_wdir, "site_summary.csv"))

# read in species metadata
species_metadata <- read_csv(here(remote_wdir, "species_metadata.csv"))

## =======================================================================
## ---------------------------- Wrangle Data -----------------------------
## =======================================================================

# ...... update rwi data ......
rwi_clean <- rwi_long %>% 
  
  # remove cols not being used to speed up join 
  select(-c(rwl, rwi_ar, rwi_nb)) %>% 
  
  # filter out predictive sample years (2019-2065)
  filter(year < 2019)


# ...... update species metadata ......
species_metadata <- species_metadata %>% 
  rename(sp_code = species_id)


# ..... update site smry data .....
site_smry <- site_smry %>%
  select(collection_id, sp_id) %>% 
  mutate(location_id = collection_id) %>%
  
  # change sp_id to lower case
  mutate(sp_code = tolower(sp_id)) %>% 
  
  # remove sp_id column
  select(-sp_id)

## =======================================================================
## ----------------------------- Join Data -------------------------------
## =======================================================================

# ..... join rwi, metadata, and site summary .....
rwi_species <- rwi_clean %>% 
  
  # join rwi to site_summary by collection id
  left_join(y = site_smry, by = 'collection_id') %>% 
  
  #join with species metadata
  left_join(y = species_metadata, by = "sp_code")

## =======================================================================
## --------------------------- Chart: Top 20 Spp -------------------------
## =======================================================================

# ..... create df listing individual trees .....
indv_tree_spp <- rwi_species %>% 
  group_by(tree) %>% 
  summarise(sp_code = first(sp_code), 
            spp = first(spp),
            common_name = first(common_name))

# ..... count occurence of each species .....
sp_count <- plyr::count(indv_tree_spp$spp) %>% 
  rename(spp = x)

# plot top 20 collection ID counts
sp_count %>% 
  drop_na(spp) %>% 
  slice_max(order_by = freq, n = 20) %>% 
  ggplot() +
  geom_col(aes(x = fct_reorder(spp, freq), y = freq), fill = "palegreen4") +
  labs(x = "Genus Species",
       y = "Number of Trees Sampled",
       title = "20 Most Prominent Species in Tree Core Dataset") +
  geom_text(aes(x = fct_reorder(spp, freq), y = freq,
                label = freq), hjust = 1.2, color = "white", size = 3.5) +
  coord_flip() +
  theme_minimal() #+
# annotate("rect", 
#          xmin = 10.5, xmax = 6.5, 
#          ymin = -1, ymax = 731, 
#          fill = NA, color = "maroon",
#          linewidth = 1.5)

