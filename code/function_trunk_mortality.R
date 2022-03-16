## Function to use in the simulations of thermal mortality

## INPUT
### tl: tolerance landscape of the simulated data of larval time of death
### sensdata: microclimatic data
### thres: threshold, the temperature at which we want to truncate thermal profiles

## OUTPUT
### final_df: a data.frame with thermal mortality for this simulated data at minute resolution

trunk_mortality <- function(tl, sensdata, thres){
  # t1 <- proc.time()
  
  sensdata <- sensdata %>% 
    map(mutate,
        day = floor((hour-1)/24) + 1,
        ta.trunk = if_else(ta.min >= thres, thres, ta.min)) %>%
    map(~ split(., .$day))
  
  daily.mort.pr.trunk <- sensdata %>%
    map(~ future_map(.x,
                     ~ dynamic.landscape.mod(ta = .x$ta.trunk,
                                             tolerance.landscape = tl[["PR"]],
                                             plot = F)))
  
  # t2 <- proc.time()- t1
  # print(paste("PR curves for the trunkated series of sensor data at", thres, "ºC required", round(t2[3]/60), "min"))
  # t3 <- proc.time()
  
  daily.mort.pn.trunk <- sensdata %>%
    map(~ future_map(.x,
                     ~ dynamic.landscape.mod(ta = .x$ta.trunk,
                                             tolerance.landscape = tl[["PN"]],
                                             plot = F)))
  
  
  # t4 <- proc.time()- t3
  # print(paste("PN curves for the trunkated series of sensor data at", thres, "ºC required", round(t4[3]/60), "min"))
  
  
  final_df <- list(PN = daily.mort.pn.trunk,
                   PR = daily.mort.pr.trunk) %>%
    map(map,
        bind_rows,
        .id = "day") %>%
    map(bind_rows,
        .id = "sensor") %>%
    bind_rows(.id = "sp") %>%
    mutate(sp_comp = if_else(sp == "PN", "P. napi", "P. rapae"))
  
  
  return(final_df)
  

}