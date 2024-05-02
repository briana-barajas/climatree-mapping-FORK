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