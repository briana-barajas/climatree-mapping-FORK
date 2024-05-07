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