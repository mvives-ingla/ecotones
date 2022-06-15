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

## Testing more factors than temperature
mods <- realtdtdata %>% 
  nest(tdt_data = -sp) %>% 
  mutate(model = map(tdt_data,
                     ~ lmer(log10(aprox_minute_dead) ~ 
                              SENSOR_mean_temp + initial.weight + Site + (1|Site:family),
                          data = .)),
         summ = map(model, summary),
         eff = map(model, car::Anova,type = "III"))


mods$summ %>% 
  set_names(nm = mods2$sp)

mods$eff %>% 
  set_names(nm = mods2$sp)


realtdtdata %>% 
  split(.$sp) %>% 
  map(~ ggplot(data = .,
               aes(x = initial.weight, y = log10(aprox_minute_dead),
             color = as.factor(treatment))) +
        geom_point() +
        geom_smooth(method = "lm") +
    labs(title = .$sp) +
  theme_bw())




