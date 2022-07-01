#! /usr/bin/env Rscript

######################################################################################
## Simulations of thermal death time data & prediction of derived thermal mortality ##
######################################################################################

# Loop with 100 simulations for different thermal thresholds

# output:
# * mort_min_tdtsim_sensors_full: predicted thermal mortality at minute resolution
# * simdata_full: data of larval death time of all the simulations
# * tdtsim_full: tolerance landscape (TDT curves) of all the simulations
# * sensors_minute_spline_: thermal profiles form microclimatic data at minute resolution



# Pacakges ----------------------------------------------------------------
library(data.table)
library(tidyverse)
library(readxl)
library(lubridate)
library(zoo)
library(colorspace)
library(patchwork)
library(future)
library(parallel)
library(furrr)



# Parallellization setup --------------------------------------------------
plan(multicore,  workers = availableCores("multicore")-40)



# Functions ---------------------------------------------------------------
source("code/Thermal_landscape_functions_mod.R")
source("code/function_trunk_mortality.R")



# Data --------------------------------------------------------------------

tdt <- read.table("data/TDT_experiment.txt", header = T, sep = ",", dec = ",")

sensor <- read.csv ("data/sensors_hourly.csv")


all.sens <- sensor %>%
  mutate(winter_jday = if_else(winter_year == 2016, winter_jday - 1, as.double(winter_jday))) %>% 
  #To guarantee that day correspondence of 2016 is the same as other years
  group_by(sensor, winter_jday, winter_hour) %>%
  summarise(ta.h = mean(TEMP, na.rm = T)) %>%
  #hourly means across years (as Rezende)
  split(.$sensor) %>%
  map(.,
      ~ spline(.$ta.h, n = (nrow(.)-1)*60)) %>%
  #Last hour (23 to 24h) of the last days can't be added because we don't know the last temperature (00h)
  map(bind_cols) %>%
  map(rename,
      hour = x,
      ta.min = y)

juliandays <- sensor %>%
  mutate(winter_jday = if_else(winter_year == 2016, winter_jday - 1, as.double(winter_jday))) %>% #Pq el 2016 tingui la mateixa correspondència de dies que els altres anys
  group_by(sensor, winter_jday, winter_hour) %>%
  summarise(ta.h = mean(TEMP, na.rm = T)) %>%
  group_by(sensor) %>%
  distinct(winter_jday) %>%
  mutate(day = min_rank(winter_jday),
         duration = n_distinct(winter_jday))



# Simulation of thermal death time and calculation of thermal mort --------

## TDT curves from original data
tdt.split <- tdt %>% 
  filter(!is.na(aprox_minute_dead),
         !is.na(SENSOR_mean_temp)) %>% 
  split(.$sp)


tdt_curve <- tdt.split %>% 
  map(~ tdt.curve(ta = .$SENSOR_mean_temp, time = .$aprox_minute_dead))



## Simulations & estimates of daily mortality

thres <- c(100, 45, 42.5, 40, 37.5, 35)
sims <- 100

pb <- txtProgressBar(min = 0, max = sims*length(thres), initial = 0, char = "*",  style = 3)

for (i in seq_along(thres)) {
  for (j in 1:sims) {
    
    ## Tolerance landscape of simulated data
    simdata <- tdt_curve %>% 
      map(pluck,
          "model") %>% 
      map(~ data.frame(treatment = rep(c(40, 42, 44), each = 35),
                       inter = unname(coefficients(.)[1], force = T),
                       slope = unname(coefficients(.)[2], force = T),
                       sigma = sigma(.))) %>% 
      map(mutate,
          error = rnorm(n = n(), sd = sigma),
          log10time = inter + slope*treatment + error,
          time = 10^log10time)
    
    
    
    tl.sim <- simdata %>% 
      map(~ tolerance.landscape(ta = .$treatment, time = .$time))
    
    
    ## Estimates of daily mortality
    
    subtibble <- trunk_mortality(tl = tl.sim, sensdata = all.sens, thres = thres[i]) %>% 
      mutate(alive = if_else(is.na(alive), 0, alive),
             mort = 100-alive,
             day = as.numeric(day)) %>% 
      left_join(juliandays,
                by = c("sensor" = "sensor", "day" = "day")) %>% 
      group_by(sensor, winter_jday, sp_comp) %>% 
      summarise(mort = max(mort, na.rm = T),
                tmax = max(ta, na.rm = T),
                tmin = min(ta, na.rm = T),
                tmean = mean(ta, na.rm = T),
                surv = min(alive, na.rm = T),
                thres = thres[i],
                simulation = j)
    
    if(i == 1 & j == 1){
      fwrite(subtibble,
             file = paste0("data/sim_mortalities/mort_min_tdtsim_sensors_full",
                           Sys.Date(),
                           ".csv"))
      
      simdata %>% 
        bind_rows(.id = "sp") %>% 
        select(sp, treatment, time) %>% 
        mutate(thres = thres[i],
               simulation = j) %>% 
        fwrite(file = paste0("data/sim_mortalities/simdata_full",
                             Sys.Date(),
                             ".csv"))
      
      tl.sim %>%
        map(~ keep(., names(.) %in% c("ctmax", "z"))) %>% 
        map(bind_rows) %>% 
        bind_rows(.id = "sp") %>% 
        mutate(simulation = j,
               thres = thres[i]) %>% 
        fwrite(file = paste0("data/sim_mortalities/tdtsim_full",
                             Sys.Date(),
                             ".csv"))
      
    } else {
      fwrite(subtibble,
             file = paste0("data/sim_mortalities/mort_min_tdtsim_sensors_full",
                           Sys.Date(),
                           ".csv"),
             append = T)
      
      simdata %>% 
        bind_rows(.id = "sp") %>% 
        select(sp, treatment, time) %>% 
        mutate(thres = thres[i],
               simulation = j) %>% 
        fwrite(file = paste0("data/sim_mortalities/simdata_full",
                             Sys.Date(),
                             ".csv"),
               append = T)
      
      tl.sim %>%
        map(~ keep(., names(.) %in% c("ctmax", "z"))) %>% 
        map(bind_rows) %>% 
        bind_rows(.id = "sp") %>% 
        mutate(simulation = j,
               thres = thres[i]) %>% 
        fwrite(file = paste0("data/sim_mortalities/tdtsim_full",
                             Sys.Date(),
                             ".csv"),
               append = T)
    }
    
    setTxtProgressBar(pb, sims*(i-1)+j)
  }
}



all.sens %>% 
  bind_rows(.id = "sensor") %>% 
  write_rds(file = paste0("data/sim_mortalities/sensors_minute_spline_",
                          Sys.Date(),
                          ".RDS"))




