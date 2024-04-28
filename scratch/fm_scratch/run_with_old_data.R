#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 5/1/20
# Purpose: 1) Characterize climate niche for different species. 
#          2) Standardize annual site data, baseline raster, and CMIP predictions for each species
#
# Input files:
#   merged_ranges.shp:
#   HistoricCWD_AETGrids_Annual.Rdat: Rasters describing historic CWD and AET
#     Generated using historic_cwdraster.R
#   monthlycrubaseline_tas:
#   cmip5_cwdaet_start.Rdat:
#   cmip5_cwdaet_end.Rdat:
#   essentialcwd_data.csv:
#   site_summary.csv:
# 
# Output files:
#   clim_niche.csv: Tabulation of each species' historic climate niche. Parameters used for standardization.
#   sp_clim_predictions.gz: Dataset describing species' standardized historic climate and their predicted climate under CMIP5 across species' full range
#   site_ave_clim.gz: Tables describing each site's species-standardized average historic climate
#   site_an_clim.gz: Tables describing each site's annual species-standardized weather
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyverse)
library(dbplyr)
library(RSQLite)
library(ggplot2)
library(sf)
library(rgeos)
library(stringr)
library(raster)
library(readr)
library(tmap)
library(tictoc)
select <- dplyr::select


library(furrr)
n_cores <- 8
future::plan(multisession, workers = n_cores)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load data --------------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
data_dir <- "~/../../capstone/climatree/raw_data/"
output_dir <- "~/../../capstone/climatree/output/1-process-raw-data/"

# 1. Historic climate raster
clim_file <- paste0(data_dir, 'HistoricCWD_AETGrids_Annual.Rdat')
load(clim_file)
cwd_historic <- mean(cwd_historic)
aet_historic <- mean(aet_historic)
pet_historic <- aet_historic + cwd_historic
names(cwd_historic) = "cwd"
names(pet_historic) = "pet"

# 2. Data on historic baseline temp and precip
temps_historic <- raster(paste0(data_dir, "monthlycrubaseline_tas"))
names(temps_historic) = "temp"
temps_historic <- resample(temps_historic, cwd_historic)
clim_historic <- raster::brick(list(cwd_historic, pet_historic, temps_historic))

# 3. Site-specific historic climate data
site_clim_csv <- paste0(data_dir, 'essentialcwd_data_old.csv')
site_clim_df <- read_csv(site_clim_csv)
site_clim_df <- site_clim_df %>% 
  mutate("site_id" = as.character(site)) %>% 
  rename(location_id = site_id)#,
        # precip = ppt) 

# 4. Load species information for sites
site_smry <- read_csv(paste0(data_dir, 'site_summary.csv'))
site_smry <- site_smry %>%
  select(collection_id, sp_id) %>% 
  mutate(location_id = collection_id) %>% 
  mutate(sp_code = tolower(sp_id)) %>% 
  select(-sp_id)


## NOTE: FIA data not included in replication data repository
# site_smry_fia <- read_csv(paste0(wdir, 'out/dendro/site_summary_fia.csv'))
# site_smry_fia <- site_smry_fia %>% 
#   select(collection_id, location_id = plot_cn, sp_id = species_id) %>% 
#   mutate(sp_code = tolower(sp_id)) %>% 
#   select(-sp_id)
# site_smry <- rbind(site_smry, site_smry_fia)


# 5. Species range maps
range_file <- paste0(data_dir, 'merged_ranges_dissolve.shp')
range_sf <- st_read(range_file)

# 6. Climate projections from CMIP5
cmip_end <- load(paste0(data_dir, 'cmip5_cwdaet_end.Rdat'))
pet_cmip_end <- aet_raster + cwd_raster
cwd_cmip_end <- cwd_raster
names(cwd_cmip_end) <- NULL # Resetting this due to strange names in file from CMIP processing
rm(cwd_raster)
rm(aet_raster)

cmip_start <- load(paste0(data_dir, 'cmip5_cwdaet_start.Rdat'))
pet_cmip_start <- aet_raster + cwd_raster
cwd_cmip_start <- cwd_raster
names(cwd_cmip_start) <- NULL # Resetting this due to strange names in file from CMIP processing
rm(cwd_raster)
rm(aet_raster)




# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Visually inspect data -----------------------------------------------
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# tmap_mode("view")
# tm_shape(cwd_cmip_end) +
#   tm_raster() +
#   tm_facets(as.layers = TRUE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summarize species niches -----------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Pull and organize climate distribution for species
pull_clim <- function(spp_code){
  print(spp_code)
  # Pull relevant range map
  sp_range <- range_sf %>%
    filter(sp_code == spp_code) %>% 
    rasterize(cwd_historic, getCover=TRUE)
  sp_range[sp_range==0] <- NA
  
  # Pull cwd and aet values
  cwd_vals <- cwd_historic %>% 
    mask(sp_range) %>% 
    as.data.frame(xy = TRUE) %>% 
    drop_na()
  
  pet_vals <- pet_historic %>% 
    mask(sp_range) %>% 
    as.data.frame(xy = TRUE) %>% 
    drop_na()
  
  temp_vals <- temps_historic %>% 
    mask(sp_range) %>% 
    as.data.frame(xy = TRUE) %>% 
    drop_na()
  
  # Combine into tibble
  clim_vals <- cwd_vals %>% 
    left_join(pet_vals, by = c("x", "y")) %>% 
    left_join(temp_vals, by = c("x", "y"))
  
  return(clim_vals)
}


species_list <- range_sf %>%
  pull(sp_code) %>% 
  unique() %>% 
  enframe(name = NULL) %>% 
  select(sp_code = value) %>% 
  arrange(sp_code) %>% 
  drop_na()

clim_df <- species_list %>% 
  mutate(clim_vals = future_map(sp_code, 
                                .f = pull_clim,
                                .options = furrr_options(packages = c( "dplyr", "raster", "sf")),
                                .progress = TRUE))


## Summarize mean and sd of each species' climate
niche_df <- clim_df %>% 
  unnest(clim_vals) %>% 
  group_by(sp_code) %>% 
  summarize(pet_mean = mean(pet),
            pet_sd = sd(pet),
            cwd_mean = mean(cwd),
            cwd_sd = sd(cwd),
            temp_mean = mean(temp),
            temp_sd = sd(temp))


## Export species niche description
write.csv(niche_df, paste0(output_dir, "clim_niche_old.csv"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Standardize historic climate -------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sp_standardize <- function(val, sp_mean, sp_sd){
  std_val <- (val - sp_mean) / sp_sd
  return(std_val)
}

sp_std_historic_df <- function(hist_clim_vals, pet_mean, pet_sd, cwd_mean, cwd_sd, temp_mean, temp_sd){
  hist_clim_vals <- hist_clim_vals %>% 
    mutate_at(vars(starts_with("cwd")), 
              ~sp_standardize(.x, cwd_mean, cwd_sd)) %>% 
    mutate_at(vars(starts_with("pet")), 
              ~sp_standardize(.x, pet_mean, pet_sd)) %>% 
    mutate_at(vars(starts_with("temp")), 
              ~sp_standardize(.x, temp_mean, temp_sd))
  return(hist_clim_vals)
}


sp_std_future_df <- function(cmip_df, hist_clim_vals, pet_mean, pet_sd, cwd_mean, cwd_sd, temp_mean, temp_sd){
  valid_locations <- hist_clim_vals %>% select(x,y)
  cmip_df <- valid_locations %>% 
    left_join(cmip_df, by = c("x", "y"))
  cmip_df <- cmip_df %>% 
    mutate_at(vars(starts_with("cwd")), 
              ~sp_standardize(.x, cwd_mean, cwd_sd)) %>% 
    mutate_at(vars(starts_with("pet")), 
              ~sp_standardize(.x, pet_mean, pet_sd)) %>% 
    mutate_at(vars(starts_with("temp")), 
              ~sp_standardize(.x, temp_mean, temp_sd))
  return(cmip_df)
}


clim_df <- clim_df %>% 
  left_join(niche_df, by = "sp_code")

clim_df <- clim_df %>% 
  mutate(clim_historic_sp = future_pmap(list(hist_clim_vals = clim_vals,
                                             pet_mean = pet_mean,
                                             pet_sd = pet_sd,
                                             cwd_mean = cwd_mean,
                                             cwd_sd = cwd_sd,
                                             temp_mean = temp_mean,
                                             temp_sd = temp_sd),
                                        .f = sp_std_historic_df,
                                        .options = furrr_options(packages = c( "dplyr"))))
# NOTE: May no longer need this dataframe???


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Pull CMIP projections -----------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pet_end_df <- pet_cmip_end %>% 
  as.data.frame(xy = TRUE) %>%
  drop_na() %>% 
  as_tibble() %>%
  rename_with(~stringr::str_replace(., "layer.", "pet_cmip_end"), 
              starts_with('layer.'))

cwd_end_df <- cwd_cmip_end %>% 
  as.data.frame(xy = TRUE) %>%
  drop_na() %>% 
  as_tibble() %>%
  rename_with(~stringr::str_replace(., "layer.", "cwd_cmip_end"), 
              starts_with('layer.'))


pet_start_df <- pet_cmip_start %>% 
  as.data.frame(xy = TRUE) %>%
  drop_na() %>% 
  as_tibble() %>%
  rename_with(~stringr::str_replace(., "layer.", "pet_cmip_start"), 
              starts_with('layer.'))

cwd_start_df <- cwd_cmip_start %>% 
  as.data.frame(xy = TRUE) %>%
  drop_na() %>% 
  as_tibble() %>%
  rename_with(~stringr::str_replace(., "layer.", "cwd_cmip_start"), 
              starts_with('layer.'))


# ## Illustrate forawrd/backward conversion between df and raster
# crs_template <- crs(cwd_future)
# raster_template <- cwd_df %>% select(x,y)
# cwd_df <- cwd_df %>% 
#   drop_na()
# cwd_df2 <- raster_template %>% 
#   left_join(cwd_df, by = c("x", "y"))
# cwd_rast2 <- rasterFromXYZ(cwd_df2, crs = crs_template)

## Combine PET and CWD projections
cmip_df <- cwd_end_df %>% 
  full_join(pet_end_df, by = c("x", "y")) %>% 
  full_join(cwd_start_df, by = c("x", "y")) %>% 
  full_join(pet_start_df, by = c("x", "y"))

## Nest CMIP data
cmip_df <- cmip_df %>%
  mutate(idx = 1) %>% 
  group_by(idx) %>% 
  nest() %>% 
  ungroup() %>% 
  select(data)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Summarize cmip climate for each species ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Cross species list with nested cmip data
sp_cmip_clim <- clim_df %>% 
  mutate(cmip_df = cmip_df$data)



sp_cmip_clim <- sp_cmip_clim %>% 
  mutate(clim_cmip_sp = future_pmap(list(cmip_df = cmip_df,
                                         hist_clim_vals = clim_vals,
                                         pet_mean = pet_mean,
                                         pet_sd = pet_sd,
                                         cwd_mean = cwd_mean,
                                         cwd_sd = cwd_sd), 
                                    .f = sp_std_future_df,
                                    .options = furrr_options(packages = c( "dplyr")))) %>% 
  select(-cmip_df)


# ## Check final result as raster
# species = "acsh"
# test_clim <- (sp_cmip_clim %>% filter(sp_code == species) %>% pull(clim_cmip_sp))[[1]]
# crs_template <- crs(cwd_cmip_end)
# raster_template <- cwd_cmip_end %>% as.data.frame(xy = TRUE) %>% select(x,y)
# test_clim <- raster_template %>%
#   left_join(test_clim, by = c("x", "y"))
# test_clim <- rasterFromXYZ(test_clim, crs = crs_template)
# range <- range_sf %>% filter(sp_code == species)
# tmap_mode("view")
# 
# tm_shape(test_clim$cwd_cmip_end1) +
#   tm_raster(palette = "-RdYlGn") +
#   tm_facets(as.layers = TRUE) +
#   tm_shape(range) + 
#   tm_fill(col = "lightblue")
# 
# tm_shape(test_clim$cwd_cmip_end1) +
#   tm_raster(palette = "-RdYlGn") +
#   tm_facets(as.layers = TRUE) +
#   tm_shape(test_clim$cwd_cmip_start1) +
#   tm_raster(palette = "-RdYlGn") +
#   tm_facets(as.layers = TRUE) +
#   tm_shape(range) + 
#   tm_fill(col = "lightblue")


## Export predictions
write_rds(sp_cmip_clim, paste0(output_dir, "sp_clim_predictions_old.", compress = "gz"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Apply species standardization to site-level data -----------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# niche_df <- niche_df %>% 
#   select(sp_code = sp_code, sp_pet_mean = pet_mean, sp_pet_sd = pet_sd, sp_cwd_mean = cwd_mean, sp_cwd_sd = cwd_sd)

# ===============================================================
# Left off here, new data does not have the same column names
# ===============================================================

# Calculate site-level annual climate
site_clim_df = site_clim_df %>%
  group_by(location_id, year) %>%
  summarise(
    aet.an = sum(aet),
    cwd.an = sum(cwd),
    pet.an = sum((aet+cwd)),
    temp.an = mean(tmean),
    .groups = "drop")

### Calculate site-level, average, historic, relative climate (for second stage)
## TODO: Note - dropping CANA323 because it has null climate data for a few months each year. might want to dig into this
site_clim_df <- site_clim_df %>% 
  filter(location_id != "CANA323")

ave_site_clim_df <- site_clim_df %>% 
  filter(year < 1980) %>% 
  group_by(location_id) %>% 
  summarise(cwd.ave = mean(cwd.an),
            pet.ave = mean(pet.an),
            cwd.sd = sd(cwd.an),
            pet.sd = sd(pet.an),
            temp.ave = mean(temp.an),
            temp.sd = sd(temp.an)) %>% 
  ungroup()

spstd_site_clim_df <- site_smry %>% 
  left_join(ave_site_clim_df, by = "location_id") %>% 
  group_by(sp_code) %>% 
  nest(data = c(collection_id, 
                cwd.ave, 
                # pet.ave, 
                temp.ave)) %>% 
  left_join(niche_df, by = ("sp_code")) %>%
  drop_na() # Dropping some species due to NA niche data

spstd_site_clim_df <- spstd_site_clim_df %>% 
  mutate(site_clim = future_pmap(list(hist_clim_vals = data,
                                      pet_mean = pet_mean,
                                      pet_sd = pet_sd,
                                      cwd_mean = cwd_mean,
                                      cwd_sd = cwd_sd,
                                      temp_mean = temp_mean,
                                      temp_sd = temp_sd),
                                 .f = sp_std_historic_df,
                                 .options = furrr_options(packages = c( "dplyr"))))

spstd_site_clim_df <- spstd_site_clim_df %>% 
  unnest(site_clim) %>% 
  rename(cwd.spstd = cwd.ave, 
         # pet.spstd = pet.ave, 
         temp.spstd = temp.ave) %>% 
  mutate(cwd.sd = cwd.sd / cwd_sd,
         # pet.sd = pet.sd / pet_sd,
         temp.sd = temp.sd / temp_sd) %>% 
  ungroup() %>% 
  select(collection_id, location_id, cwd.spstd, 
         # pet.spstd, 
         temp.spstd, cwd.sd, 
         # pet.sd, 
         temp.sd)

spstd_site_clim_df <- spstd_site_clim_df %>% 
  left_join(ave_site_clim_df %>% select(location_id, 
                                        cwd.ave, 
                                        # pet.ave, 
                                        temp.ave), by = "location_id")

spstd_site_clim_df <- spstd_site_clim_df %>% 
  select(-location_id)

write_rds(spstd_site_clim_df, 
          paste0(output_dir, "site_ave_clim_old.", compress = "gz"))





### Calculate site-level, annual, historic, relative climate (for first stage) 
an_site_clim_df <- site_smry %>% 
  left_join(site_clim_df, by = "location_id") %>% 
  group_by(sp_code) %>% 
  nest() %>% 
  left_join(niche_df, by = "sp_code") %>% 
  drop_na()

an_site_clim_df <- an_site_clim_df %>% 
  mutate(site_clim = future_pmap(list(hist_clim_vals = data,
                                      pet_mean = pet_mean,
                                      pet_sd = pet_sd,
                                      cwd_mean = cwd_mean,
                                      cwd_sd = cwd_sd,
                                      temp_mean = temp_mean,
                                      temp_sd = temp_sd),
                                 .f = sp_std_historic_df,
                                 .options = furrr_options(packages = c( "dplyr"))))

an_site_clim_df <- an_site_clim_df %>% 
  unnest(site_clim) %>% 
  rename(cwd.an.spstd = cwd.an, 
         # pet.an.spstd = pet.an, 
         temp.an.spstd = temp.an) %>% 
  ungroup() %>% 
  select(
    # -aet.an, 
    -pet_mean, 
    -pet_sd, -cwd_mean, -cwd_sd, -temp_mean, -temp_sd, -data, -sp_code)

an_site_clim_df <- an_site_clim_df %>%
  select(-location_id)

write_rds(an_site_clim_df, 
          paste0(output_dir, "site_an_clim_old.", compress = "gz"))


# ## Exploring source of dropped sites - seems to be entirely driven by sites for species with no range maps
# an_site_clim_df %>% pull(collection_id) %>% unique() %>% length()
# site_clim_df %>% pull(collection_id) %>% unique() %>% length()
# # clim_sites <- clim_df %>% pull(collection_id) %>% unique()
# test_sites <- test %>% pull(collection_id) %>% unique()
# an_site_clim_df %>% unnest(cols = c(data)) %>% pull(collection_id) %>% unique() %>% length()
# an_site_clim_df %>% unnest(cols = c(data)) %>% drop_na() %>% pull(collection_id) %>% unique() %>% length()








4A 
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
4A










#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Robert Heilmayr, Frances Moore, Joan Dudney
# Project: Treeconomics
# Date: 5/1/20
# Purpose: Run plot-level regressions of RWI sensitivity to annual weather variability
#
# Input files:
# - site_an_clim.gz: File detailing site-level weather history. Generated using "3b. Species_niche.R"
# - rwi_long.csv: Data containing processed RWI data. Generated using "1b. Parse ITRDB.R"
# - site_summary.csv: Summary data about each site. Generated using "1b. Parse ITRDB.R"
#
# Output files:
# - example_sites.csv: Dendrochronologies for two example sites. Used for methods summary figure.
# - site_pet_cwd_std.csv: Table of first stage regression parameters for baseline specification.
# - site_pet_cwd_std_nb.csv: Table of first stage regression parameters for robustness model using NB desplining.
# - site_pet_cwd_std_ar.csv: Table of first stage regression parameters for robustness model using AR desplining.
# - site_temp_cwd_std.csv:Table of first stage regression parameters for robustness model using temperature in place of PET.
#
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyr)
library(tidyverse)
# library(tidylog)
library(dbplyr)
library(broom.mixed)
library(broom)
library(purrr)
library(fixest)
library(dtplyr)
library(furrr)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import and integrate data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
data_dir <- "~/../../capstone/climatree/raw_data/"
output_dir <- "~/../../capstone/climatree/output/1-process-raw-data/"


# 1. Dendrochronologies
#dendro_dir <- paste0(wdir, "1_input_processed/dendro/")
dendro_df <- read_csv(paste0(output_dir, "rwi_long_subset.csv"))
dendro_df <- dendro_df %>% 
  select(-core_id)

## Combine multiple cores from the same tree
dendro_df <- dendro_df %>% 
  lazy_dt() %>% 
  group_by(collection_id, tree, year) %>% 
  summarise(rwi = mean(rwi),
            rwl = mean(rwl),
            rwi_ar = mean(rwi_ar),
            rwi_nb = mean(rwi_nb),
            .groups = "drop") %>% 
  as_tibble()

# 2. Historic site-level climate
an_site_clim <- read_rds(paste0(output_dir, "site_an_clim_old.gz"))
dendro_df <- dendro_df %>% 
  left_join(an_site_clim, by = c("collection_id", "year"))


# 3. Site information
site_smry <- read_csv(paste0(data_dir, 'site_summary.csv'))
site_smry <- site_smry %>% 
  select(collection_id, sp_id) %>% 
  mutate(species_id = tolower(sp_id)) %>% 
  select(-sp_id)

dendro_df <- dendro_df %>% 
  left_join(site_smry, by = 'collection_id')


# 4. Drop data from species without range maps and resulting climatic niche data
niche_df <- read.csv(paste0(output_dir, "clim_niche.csv")) %>%
  select(-X)
niche_species <- niche_df %>% pull(sp_code) %>% unique()
dendro_species <- dendro_df %>% pull(species_id) %>% unique()
dendro_df <- dendro_df %>% 
  filter(species_id %in% niche_species)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Export example sites for presentations  ------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ex_sites <- c("CO559", "CA585")
# dendro_df %>% 
#   filter(collection_id %in% ex_sites) %>% 
#   write.csv(paste0(output_dir, "example_sites.csv"))



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define regression model  -------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs_mod <- function(site_data, outcome = "rwi", energy_var = "pet.an", mod_type = "lm"){
  failed <- F
  reg_error <- NA
  nobs <- NA
  ntrees <- site_data %>% select(tree) %>%  n_distinct()
  no_cwd_var <- (site_data %>% select(cwd.an) %>% n_distinct() == 1)
  no_pet_var <- (site_data %>% select(energy_var) %>% n_distinct() == 1)
  
  if (no_cwd_var | no_pet_var) {
    message(paste0("Site has no variation in cwd.an or ", energy_var))
    failed <- T
  } else{
    # Try to run felm. Typically fails if missing cwd / pet data 
    tryCatch(
      expr = {
        formula <- as.formula(paste0(outcome, " ~ ", energy_var, " + cwd.an"))
        if (mod_type == "lm"){
          mod <- lm(formula, data = site_data)
        }
        if (mod_type == "lme"){
          mod <- nlme::lme(formula,
                           data=site_data, method="REML",
                           random = ~ 1 | tree,
                           correlation = nlme::corAR1(form=~year|tree))
        }
        
        mod_sum <- summary(mod)
        mod_vcov <- vcov(mod)
        # cov <- list(int_cwd = mod_vcov[1, 2], 
        #             int_pet = mod_vcov[1, 3], 
        #             pet_cwd = mod_vcov[2, 3])
        nobs <- nobs(mod)
        mod <- tidy(mod) %>%
          mutate(term = term %>% str_replace("\\(Intercept\\)", "intercept")) %>% 
          filter(term %in% c('intercept', 'cwd.an', energy_var)) %>% 
          pivot_wider(names_from = "term", values_from = c("estimate", "std.error", "statistic", "p.value"))
        # mod <- mod %>% 
        #   rename_all(funs(stringr::str_replace_all(., energy_var, 'energy.an')))
        mod$cov_int_cwd = mod_vcov[c("(Intercept)"), c("cwd.an")]
        cov_var_name <- paste0("cov_int_", energy_var %>% str_replace(".an", ""))
        mod[[cov_var_name]] = mod_vcov[c("(Intercept)"), c(energy_var)]
        cov_var_name <- paste0("cov_cwd_", energy_var %>% str_replace(".an", ""))
        mod[[cov_var_name]] = mod_vcov[c("cwd.an"), c(energy_var)]
        mod$r2 = mod_sum$r.squared
      },
      error = function(e){ 
        message("Returned regression error")
        reg_error <<- e[1]
        failed <<- T
      }
    )    
  }
  if (failed){
    return(NA)
  }
  return(tibble(mod = list(mod), nobs = nobs, ntrees = ntrees, error = reg_error))
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run site-level regressions --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
site_df <- dendro_df %>% 
  # drop_na() %>% 
  rename(cwd.an = cwd.an.spstd,
         #pet.an = pet.an.spstd,
         temp.an = temp.an.spstd) %>% 
  group_by(collection_id) %>%
  add_tally(name = 'nobs') %>% 
  # filter(nobs>10) %>% 
  nest()


fs_mod_bl <- partial(fs_mod, outcome = "rwi", energy_var = "pet.an", mod_type = "lm")
# fs_mod_nb <- partial(fs_mod, outcome = "rwi_nb", energy_var = "pet.an", mod_type = "lm")
# fs_mod_ar <- partial(fs_mod, outcome = "rwi_ar", energy_var = "pet.an", mod_type = "lm")
# fs_mod_temp <- partial(fs_mod, outcome = "rwi", energy_var = "temp.an", mod_type = "lm")
# fs_mod_re <- partial(fs_mod, outcome = "rwi", energy_var = "pet.an", mod_type = "lme")

site_df <- site_df %>% 
  mutate(fs_result = map(data, .f = fs_mod_bl))
         # fs_result_nb = map(data, .f = fs_mod_nb),
         # fs_result_ar = map(data, .f = fs_mod_ar),
         # fs_result_temp = map(data, .f = fs_mod_temp),
         # fs_result_re = map(data, .f = fs_mod_re))


data_df <- site_df %>% 
  select(collection_id,data)

fs_df <- site_df %>% 
  select(collection_id, fs_result) %>% 
  unnest(fs_result)

fs_df <- fs_df[which(!(fs_df %>% pull(mod) %>% is.na())),]
fs_df <- fs_df %>% 
  unnest(mod)

fs_df <- fs_df %>% 
  select(-error)

fs_df %>% write_csv(paste0(output_dir, 'site_pet_cwd_std_old.csv'))


# ## Repeat using results from nb detrended data
# fs_nb <- site_df %>% 
#   select(collection_id, fs_result_nb) %>% 
#   unnest(fs_result_nb)
# fs_nb <- fs_nb[which(!(fs_nb %>% pull(mod) %>% is.na())),]
# fs_nb <- fs_nb %>% 
#   unnest(mod) %>% 
#   select(-error)
# fs_nb %>% write_csv(paste0(output_dir, 'site_pet_cwd_std_nb_old.csv'))
# 
# 
# ## Repeat using results from ar detrended data
# fs_ar <- site_df %>% 
#   select(collection_id, fs_result_ar) %>% 
#   unnest(fs_result_ar)
# fs_ar <- fs_ar[which(!(fs_ar %>% pull(mod) %>% is.na())),]
# fs_ar <- fs_ar %>% 
#   unnest(mod) %>% 
#   select(-error)
# fs_ar %>% write_csv(paste0(output_dir, 'site_pet_cwd_std_ar_old.csv'))
# 
# 
# ## Repeat using results from temp model
# fs_temp <- site_df %>% 
#   select(collection_id, fs_result_temp) %>% 
#   unnest(fs_result_temp)
# fs_temp <- fs_temp[which(!(fs_temp %>% pull(mod) %>% is.na())),]
# fs_temp <- fs_temp %>% 
#   unnest(mod) %>% 
#   select(-error)
# fs_temp %>% write_csv(paste0(output_dir, 'site_temp_cwd_std_old.csv'))
# 
# 
# ## Repeat using results from re model
# fs_re <- site_df %>% 
#   select(collection_id, fs_result_re) %>% 
#   unnest(fs_result_re)
# fs_re <- fs_re[which(!(fs_re %>% pull(mod) %>% is.na())),]
# fs_re <- fs_re %>% 
#   unnest(mod) %>% 
#   select(-error)
# fs_re %>% write_csv(paste0(output_dir, 'site_pet_cwd_std_re_old.csv'))






5A
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
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
flm_df1 <- read_csv(paste0(output_dir, 'site_pet_cwd_std_old.csv'))

# 2. Historic site-level climate
ave_site_clim <- read_rds(paste0(output_dir, "site_ave_clim_old.gz"))
flm_df1 <- flm_df1 %>% 
  left_join(ave_site_clim, by = c("collection_id"))

# 3. Site information
site_df <- read_csv(paste0(data_dir, 'site_summary.csv'))
site_df <- site_df %>% 
  select(collection_id, sp_id, latitude, longitude)
site_df <- site_df %>% 
  rename(species_id = sp_id) %>% 
  mutate(species_id = str_to_lower(species_id))

# # 4. Species information
 sp_info <- read_csv(paste0(data_dir, 'species_metadata.csv'))
 sp_info <- sp_info %>% 
   select(species_id, genus, gymno_angio, family)
 site_df <- site_df %>% 
   left_join(sp_info, by = "species_id")

# Merge back into main flm_df
flm_df_old <- flm_df1 %>% 
  left_join(site_df, by = "collection_id")



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

# flm_df <- flm_df %>%
#   mutate(outlier = (estimate_cwd.an<cwd_est_bounds[1]) |
#            (estimate_cwd.an>cwd_est_bounds[2]) |
#            (estimate_pet.an<pet_est_bounds[1]) |
#            (estimate_pet.an>pet_est_bounds[2]) |
#            (cwd.spstd<cwd_spstd_bounds[1]) |
#            (cwd.spstd>cwd_spstd_bounds[2]) |
#            (pet.spstd<pet_spstd_bounds[1]) |
#            (pet.spstd>pet_spstd_bounds[2]))


flm_df <- flm_df %>%
  mutate(outlier = (estimate_cwd.an<cwd_est_bounds[1]) |
           (estimate_cwd.an>cwd_est_bounds[2]) |
           (estimate_pet.an<pet_est_bounds[1]) |
           (estimate_pet.an>pet_est_bounds[2]))

# Save out full flm_df to simplify downstream scripts and ensure consistency
flm_df %>% write.csv(paste0(output_dir, "site_pet_cwd_std_augmented.csv"))

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
save(block_list,file=paste0(output_dir,"spatial_blocks.Rdat"))
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
block_draw_df %>% 
  # select(boot_id, collection_id, cwd_coef, pet_coef, int_coef, cwd.spstd, pet.spstd) %>% 
  write_rds(paste0(output_dir, "mc_sample.gz"), compress = "gz")


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
write_rds(boot_df, paste0(output_dir, "ss_bootstrap_old.rds"))










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
output_dir <- "~/../../capstone/climatree/output/1-process-raw-data/"

# Create output directories
#out_dir <- paste0(wdir,"2_output/predictions/")
# dir.create(file.path(out_dir), showWarnings = FALSE)
#dir.create(file.path(paste0(output_dir, "sp_rwi/")), showWarnings = FALSE)
# dir.create(file.path(paste0(out_dir, "sp_hot_cells/")), showWarnings = FALSE)

# 1. Second stage model
mod_df <- read_rds(paste0(output_dir, "ss_bootstrap.rds"))
mod_df <- mod_df %>% 
  rename(iter_idx = boot_id)

# 2. Species-standardized historic and future climate
sp_clim <- read_rds(paste0(output_dir, "sp_clim_predictions.gz"))
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
    write_rds(file = paste0(output_dir, "sp_rwi/", spp_code, ".gz"), compress = "gz")
  
  toc()
  return(agg_stats)
}


mc_nests <- sp_mc %>%
  group_by(sp_code) %>%
  nest() %>% 
  drop_na()

# Generally have memory issues with 38 (LAGM), and 93 (PISY) - need to run these with two cores
# large_range_sp <- c("lagm", "pisy")
# spp_code = "abal"
# mc_data = mc_nests %>% filter(sp_code == spp_code) %>% pull(data)
# mc_data = mc_data[[1]]
# parallel = FALSE

# mc_nests_large <- mc_nests %>% 
#   filter((sp_code %in% large_range_sp)) %>% 
#   mutate(predictions = pmap(list(spp_code = sp_code,
#                                  mc_data = data,
#                                  parallel = TRUE),
#                             .f = calc_rwi_quantiles))

mc_nests_small <- mc_nests %>% 
  # filter(!(sp_code %in% large_range_sp)) %>% 
  mutate(predictions = pmap(list(spp_code = sp_code,
                                 mc_data = data,
                                 parallel = TRUE),
                            .f = calc_rwi_quantiles)) 

agg_stats <- mc_nests_small %>% 
  select(-data) %>% 
  unnest(predictions) %>% 
  write_rds(file = paste0(output_dir, "mc_agg_stats.gz"), compress = "gz")


test <- agg_stats %>% 
  group_by(iter_idx) %>% 
  summarise(rwi_pred_change = mean(rwi_pred_change))
test %>%
  pull(rwi_pred_change) %>% 
  quantile(c(0.025, 0.5, 0.975))




# # Profiling of main function
# spp_code = "juex"
# mc_data = (mc_nests %>% filter(sp_code == spp_code) %>% pull(data))[[1]]
# l = profvis(calc_rwi_quantiles(spp_code, mc_data))






# %>% 
#   select(-data) %>% 
#   unnest(predictions)

# mc_nests %>% 
#   saveRDS(file = paste0(wdir,"out/predictions/rwi_predictions.rds"))


# # # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # # Thought experiments - partialling out mechanisms    ---------------
# # # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # calc_rwi_partial_sens <- function(cmip_rast, sensitivity){ # NOTE: Should the means be calculated across full range rather than by each species?
# #   mean_fut_cwd <- cmip_rast %>% subset("cwd.spstd") %>% cellStats(stat = "mean")
# #   mean_fut_pet <- cmip_rast %>% subset("pet.spstd") %>% cellStats(stat = "mean")
# #   cwd_sens = sensitivity %>% subset("cwd_sens")
# #   pet_sens = sensitivity %>% subset("pet_sens")
# #   intercept = sensitivity %>% subset("intercept")
# #   rwi_rast <- intercept + (mean_fut_cwd * cwd_sens) + (mean_fut_pet * pet_sens)
# #   names(rwi_rast) = "rwi_psens"
# #   return(rwi_rast)
# # }
# # 
# # calc_rwi_partial_clim <- function(cmip_rast, sensitivity){ # NOTE: Should the means be calculated across full range rather than by each species?
# #   mean_cwd_sens <- sensitivity %>% subset("cwd_sens") %>% cellStats(stat = "mean")
# #   mean_pet_sens <- sensitivity %>% subset("pet_sens") %>% cellStats(stat = "mean")
# #   mean_intercept <-sensitivity %>% subset("intercept") %>% cellStats(stat = "mean") 
# #   rwi_rast <- mean_intercept + (cmip_rast$cwd.spstd * mean_cwd_sens) + (cmip_rast$pet.spstd * mean_pet_sens)
# #   names(rwi_rast) = "rwi_pclim"
# #   return(rwi_rast)
# # }
# # 
# # 
# # sp_predictions <- sp_predictions %>% 
# #   mutate(rwi_predictions_partial_sens = map2(.x = clim_future_sp, .y = sensitivity, calc_rwi_partial_sens),
# #          rwi_predictions_partial_clim = map2(.x = clim_future_sp, .y = sensitivity, calc_rwi_partial_clim))
# 
# 
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # Stack rasters into dataframe ------------------------------------
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# create_prediction_df <- function(spp_predictions){
#   sp_fut <- (spp_predictions %>% 
#                pull(clim_future_sp))[[1]]
#   names(sp_fut) <- c("cwd.fut", "pet.fut")
#   
#   sp_hist <- (spp_predictions %>% 
#                 pull(clim_historic_sp))[[1]]
#   sp_sens  <- (spp_predictions %>% 
#                  pull(sensitivity))[[1]]
#   sp_rwi  <- (spp_predictions %>% 
#                 pull(rwi_predictions))[[1]]
#   
#   sp_rwi_psens  <- (spp_predictions %>% 
#                       pull(rwi_predictions_partial_sens))[[1]]
#   
#   sp_rwi_pclim  <- (spp_predictions %>% 
#                       pull(rwi_predictions_partial_clim))[[1]]
#   
#   clim_compare <- brick(c(sp_fut, sp_hist, sp_sens, sp_rwi, sp_rwi_psens, sp_rwi_pclim))
#   clim_compare <- clim_compare %>% 
#     as.data.frame(xy = TRUE) %>% 
#     drop_na()
#   return(clim_compare)
# }
# 
# sp_prediction_test <- sp_predictions %>% 
#   group_by(sp_code, iter_idx) %>% 
#   nest() %>% 
#   mutate(pred_df = map(data, create_prediction_df)) %>% 
#   select(-data) %>% 
#   unnest(cols = pred_df) %>% 
#   mutate(cwd_change = cwd.fut - cwd.spstd,
#          pet_change = pet.fut - pet.spstd,
#          rwi_null = cwd.spstd * cwd_sens + pet.spstd * pet_sens + intercept,
#          rwi_change = rwi_pred - rwi_null,
#          rwi_change_psens = rwi_psens - rwi_null,
#          rwi_change_pclim = rwi_pclim - rwi_null)
# 
# 
# sp_predictions %>% 
#   saveRDS(file = paste0(wdir,"out/predictions/sp_predictions.rds") )





# crs_template <- crs(cwd_future)
# cwd_df <- cwd_rast %>% as.data.frame(xy = TRUE) 
# raster_template <- cwd_df %>% select(x,y)
# cwd_df <- cwd_df %>% 
#   drop_na()
# 
# cwd_df2 <- raster_template %>% 
#   left_join(cwd_df, by = c("x", "y"))
# cwd_rast2 <- rasterFromXYZ(cwd_df2, crs = crs)
# 
# tm_raster(cwd_rast2)
# data("World")
# 
# tmap_mode("view")
# tm_shape(cwd_rast) +
#   tm_raster()
# 
# tmap_mode("view")
# tm_shape(cwd_vals) +
#   tm_raster()
# 
# 
# 
# 
# 
# 
# ############ FUNCTION GRAVEYARD
# 
# calc_rwi <- function(sppp_code, cmip_id, sensitivity){
#   ## Function used to predict species' RWI rasters based on predicted 
#   ## sensitivity raster and assigned CMIP model of future climate
#   
#   sp_fut_clim <- readRDS(paste0(sp_fut_clim_dir, sppp_code, ".gz"))
#   
#   cmip_rast <- sp_fut_clim %>% 
#     filter(cmip_idx == cmip_id) %>% 
#     pull(clim_future_sp)
#   
#   cwd_sens = sensitivity %>% subset("cwd_sens")
#   pet_sens = sensitivity %>% subset("pet_sens")
#   intercept = sensitivity %>% subset("intercept")
#   rwi_rast <- intercept + (cmip_rast[[1]]$cwd.spstd * cwd_sens) + (cmip_rast[[1]]$pet.spstd * pet_sens)
#   names(rwi_rast) = "rwi_pred"
#   return(rwi_rast)
# }
# 
# quantiles <- function(x){
#   ## Defines quantiles used to summarize MC runs
#   quantile(x, c(0.025, 0.975), na.rm=TRUE)
# }
# 
# # extract_quantiles <- function(rwi_preds){
# #   ## Extracts desired quantiles for MC runs
# #   
# #   rwi_quantiles <- rwi_preds %>% 
# #     pull(rwi_predictions) %>%
# #     brick() %>% 
# #     calc(quantiles)
# #   return(rwi_quantiles)
# # }
# 
# 
# calc_mean_fut_clim <- function(sppp_code){
#   sp_fut_clim <- readRDS(paste0(sp_fut_clim_dir, sppp_code, ".gz"))
#   clim_pulls <- sp_fut_clim %>%
#     mutate(cwd.fut = pmap(list(brick = clim_future_sp, layer_name = "cwd.spstd"),
#                           .f = pull_layer),
#            pet.fut = pmap(list(brick = clim_future_sp, layer_name = "pet.spstd"),
#                           .f = pull_layer))
#   
#   cwd.fut.q <- clim_pulls %>%
#     pull(cwd.fut) %>%
#     brick()
#   cwd.fut.q <- raster::mean(cwd.fut.q)
#   names(cwd.fut.q) = "cwd.fut"
#   
#   pet.fut.q <- clim_pulls %>%
#     pull(pet.fut) %>%
#     brick()
#   pet.fut.q <- raster::mean(pet.fut.q)
#   names(pet.fut.q) = "pet.fut"
#   
#   out_brick <- brick(c(cwd.fut.q, pet.fut.q))
#   return(out_brick)
# }



