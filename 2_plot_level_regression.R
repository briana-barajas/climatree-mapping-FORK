#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Authors: Briana Barajas, Fletcher McConnell, Rosemary Juarez, Vanessa Salgado
# Project: Mapping Tree Species' Drought Sensitivity Under Climate Change
# Institution: Bren School of Environmental Science & Management - UCSB
# Date: 2024-06-07
# Purpose: Estimate site-level sensitivity to variation in climatic water deficit
#
# Input files:
# - rwi_long.csv: de-trended tree ring data from the ITRDB
# - site_summary.csv
# - site_an_clim.gz
#
# Output files:
# - site_pet_cwd_std.csv
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package imports --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyverse)
library(dbplyr)
library(broom.mixed)
library(broom)
library(purrr)
library(fixest)
library(dtplyr)
library(furrr)
select <- dplyr::select


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Import and integrate data --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Define path
data_dir <- "~/../../capstone/climatree/raw_data/"
output_dir <- "~/../../capstone/climatree/output/"


# 1. Dendrochronologies
#dendro_dir <- paste0(wdir, "1_input_processed/dendro/")
dendro_df <- read_csv(paste0(data_dir, "rwi_long.csv"))
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
an_site_clim <- read_rds(paste0(output_dir, "site_an_clim.gz"))
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
#   write.csv(paste0(wdir, "2_output/first_stage/example_sites.csv"))



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
        nobs <- nobs(mod)
        mod <- tidy(mod) %>%
          mutate(term = term %>% str_replace("\\(Intercept\\)", "intercept")) %>% 
          filter(term %in% c('intercept', 'cwd.an', energy_var)) %>% 
          pivot_wider(names_from = "term", values_from = c("estimate", "std.error", "statistic", "p.value"))
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
  rename(cwd.an = cwd.an.spstd,
         pet.an = pet.an.spstd) %>% 
  group_by(collection_id) %>%
  add_tally(name = 'nobs') %>% 
  nest()


fs_mod_bl <- partial(fs_mod, outcome = "rwi", energy_var = "pet.an", mod_type = "lm")

site_df <- site_df %>% 
  mutate(fs_result = map(data, .f = fs_mod_bl))


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

fs_df %>% write_csv(paste0(output_dir, 'site_pet_cwd_std.csv'))

