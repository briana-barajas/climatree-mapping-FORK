# Pull and organize climate distribution for species
pull_clim <- function(spp_code){
  print(spp_code)
  # Pull relevant range map
  sp_range <- range_sf %>%
    filter(sp_code == spp_code) %>% 
    rasterize(cwd_historic, getCover=TRUE)
  sp_range[sp_range==0] <- NA