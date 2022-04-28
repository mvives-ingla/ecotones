##############
# TDT models #
##############

# Packages ----------------------------------------------------------------
library (tidyverse)
library(broom)


# Functions ---------------------------------------------------------------
source("code/Thermal_landscape_functions_mod.R")


# Data --------------------------------------------------------------------
## Data from TDT experiments
realtdtdata <- read.table("data/TDT_experiment.txt", header = T,
                          sep = ",", dec = ",")


# Models ------------------------------------------------------------------

## ANCOVA: species x treatment to test whether the two species present signficantly different TDT curves
ancova <- lm(log10(aprox_minute_dead) ~ sp*SENSOR_mean_temp, data = realtdtdata) %>% 
  tidy()

## TDT curves by species
tdt <- realtdtdata %>% 
  nest(tdt_data = -sp) %>% 
  mutate(tdt = map(tdt_data,
                   ~ tdt.curve(ta = .$SENSOR_mean_temp,
                               time = .$aprox_minute_dead)))

## Testing more factor than temperature
mods <- realtdtdata %>% 
  nest(tdt_data = -sp) %>% 
  mutate(model = map(tdt_data,
                     ~ lm(log10(aprox_minute_dead) ~ 
                              SENSOR_mean_temp + initial.weight + Site*family,
                          data = .)),
         summ = map(model, summary),
         eff = map(model, anova),
         gla = map(model, glance))

mods$eff
mods$gla
