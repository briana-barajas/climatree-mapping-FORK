#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Briana Barajas, Fletcher McConnell, Vanessa Salgado, Rosemary Juarez
# Project: Mapping Global Tree Vulnerability Under Climate Change
# Date: 2024-01-28
# Purpose: 1) Isolate the 20 most common tree species in the ITRDB
#          2) Check number of sites/trees for key conserevation species
#          3) Subset rwi_long.csv to species being mapped
#
# Input files:
#   site_summary.csv
#   rwi_long.csv
#   species_metadata.csv
#
# Output files:
#   rwi_long_subset.csv (trimmed version of rwi_long with key species)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## =======================================================================
## -------------------------- Package imports ----------------------------
## =======================================================================

library(tidyverse)
library(here)

## =======================================================================
## ---------------------------- Load data --------------------------------
## =======================================================================

# set data and output directories
data_dir <- here("~","..", "..", "capstone", "climatree", "raw_data")
output_dir <- here("~","..", "..", "capstone", "climatree", "output", "1-process-raw-data")

# read in tree ring data
rwi_long <- read_csv(here(data_dir, "rwi_long.csv"))

# read in site summary table
site_smry <- read_csv(here(data_dir, "site_summary.csv"))

# read in species metadata
species_metadata <- read_csv(here(data_dir, "species_metadata.csv"))

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
# metadata contains 159 species codes
species_metadata <- species_metadata %>% 
  rename(sp_code = species_id) 


# ..... update site smry data .....
# site smry has 243 species codes
site_smry <- site_smry %>%
  
  # reduce cols to facilitate join
  select(collection_id, sp_id, n_trees) %>%

  # change codes to lowercase
  mutate(sp_code = tolower(sp_id)) %>%
  
  # remove captialized ID column
  select(-sp_id)

## =======================================================================
## ----------------------- Missing Spp Metadata --------------------------
## =======================================================================
# compare spp_codes in metadata (159) vs site summary (243)
length(unique(site_smry$sp_code)) == length(unique(species_metadata$sp_code))

# total number of individual trees in rwi_long (3100) vs site_smry (98105)
length(unique(rwi_clean$tree)) == sum(site_smry$n_trees, na.rm = TRUE)

# check for missing values
colSums(is.na(rwi_clean))
colSums(is.na(site_smry))

## =======================================================================
## ----------------------------- Join Data -------------------------------
## =======================================================================

# ..... join rwi, metadata, and site summary .....
rwi_species <- rwi_clean %>% 
  
  # join rwi to site_summary by collection id
  left_join(y = site_smry, by = 'collection_id') %>% 
  
  #join with species metadata
  left_join(y = species_metadata, by = "sp_code")

# check missing values
# colSums(is.na(rwi_species))

## =======================================================================
## --------------------------- Filter Top 20 Spp -------------------------
## =======================================================================

# ..... create df listing each trees once .....
indv_tree_spp <- rwi_species %>% 
  select(-c(core, year, gymno_angio, family)) %>% 
  
  # group by individual tree sampled
  group_by(tree) %>% 
  
  # isolate first occurrence of the tree and it's species code
  summarise(tree = first(tree),
            sp_code = first(sp_code))

# ..... count occurrence of each species .....
sp_count <- plyr::count(indv_tree_spp$sp_code) %>% 
  rename(sp_code = x) %>% 
  inner_join(species_metadata, by = "sp_code") %>% 
  slice_max(order_by = freq, n = 20) %>% 
  arrange(sp_code)

## =======================================================================
## ------------------------ Compare Filter to n_trees --------------------
## =======================================================================
# top 20 using n_trees instead of tree
top_n_trees <- site_smry %>% 
  group_by(sp_code) %>% 
  summarise(n_trees = sum(n_trees, na.rm = TRUE)) %>% 
  slice_max(order_by = n_trees, n = 20) %>% 
  arrange(sp_code)

# n_trees in original site summary data is missing, so sp_count using tree will be used
top_n_trees$sp_code
sp_count$sp_code


## =======================================================================
## --------------------------- Chart: Top 20 Spp -------------------------
## =======================================================================

# plot top 20 collection ID counts
sp_count %>% 
  ggplot() +
  geom_col(aes(x = fct_reorder(spp, freq), y = freq), fill = "palegreen4") +
  labs(x = "Genus Species",
       y = "Number of Trees Sampled",
       title = "20 Most Prominent Species in Tree Core Dataset") +
  geom_text(aes(x = fct_reorder(spp, freq), y = freq,
                label = freq), hjust = 1.2, color = "white", size = 3.5) +
  coord_flip() +
  theme_minimal()


## =======================================================================
## ---------------------------- Subset RWI Data --------------------------
## =======================================================================
# .......... create subset with key species ........
# isolate top 20 species 
top_sp_codes <- sp_count$sp_code

# add species code for species of known, high concern
top_sp_codes <- c(top_sp_codes, "pilo", "piar", "pila", "pial")

# join original data to site summary, filtering based on species list
rwi_long_subset <- rwi_long %>% 
  
  # filter out predictive sample years (2019-2065)
  filter(year < 2019) %>% 
  
  # join to site summary
  left_join(site_smry, by = "collection_id") %>% 
  
  # filter to species in list
  filter(sp_code %in% top_sp_codes) %>% 
  
  # filter back to original columns
  select(-c(n_trees, sp_code))

# .......... create subset of key species metadata ........
top_sp_codes <- species_metadata %>% 
  filter(sp_code %in% top_sp_codes) %>% 
  select(-c(source, gymno_angio, family))

# .......... write csv of rwi_long subset ........

write_csv(rwi_long_subset, here(output_dir, "rwi_long_subset.csv"))

write_csv(top_sp_codes, here(output_dir, "top_sp_codes.csv"))


