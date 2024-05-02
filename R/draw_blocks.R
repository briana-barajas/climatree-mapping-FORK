draw_blocks <- function(site_sample){
  samp <- site_sample %>% pull(samp)
  blocked_draw <- (block_list[samp] %>% unlist() %>% unname())[1:n_sites]
  return(blocked_draw)  
}