############
## Fig. 2 ##
############


# Packages ----------------------------------------------------------------
library(tidyverse)
library(readxl)
library(patchwork)
library(scales)
library(ggtext)

# Functions ---------------------------------------------------------------
source("code/boxplot_ecotones_v11.R")
source("code/result_bivar_in_subtitle.R")


# Data --------------------------------------------------------------------
## Oviposition by female
ovi.byfem <- read.csv ("data/ovi_by_female.csv")
## Temperature during ovipostion
ovi.temp <- read.csv("data/ovi_temp.csv")
## Field campaign data at the host plant level
data.plants <- read.csv("data/data_plants.csv")
## Field campaign data at the foliar level
data.leaves <- read.csv("data/data_leaves.csv")

## Microhabitat air temperature at hourly resolution (all time series)
sensor <- read.csv("data/sensors_hourly.csv")

## Daily mean and maximum temperature from weather stations (macroclimate)
meteo <- read.csv("data/meteorological_data.csv")

theme_set(theme_classic())

# A: ovipositions per microhabitat -----------------------------------------
(ovi.plot <- ovi.byfem %>% 
  mutate(microhabitat = factor(microhabitat, levels = c("C", "OC", "O"))) %>% 
  group_by(site, species, microhabitat, survey) %>% 
  summarize (ovi = sum(n_ovi, na.rm = T),
             fem = n_distinct(fem_code, na.rm = T),
             duration = mean(duration, na.rm = T)) %>% 
  group_by(species, microhabitat) %>% 
  summarize(ovi = sum(ovi, na.rm = T),
            fem = sum(fem, na.rm = T),
            duration = sum(duration, na.rm = T)) %>% 
  group_by(species) %>% 
  mutate(std_ovi = ovi/duration*600, #number of ovipositions observed in 10h (oviposition rate)
         total_std_ovi = sum (std_ovi, na.rm = T),
         #relative fraction of oviposition rate
         rel_std_ovi = std_ovi/total_std_ovi) %>% 
  ggplot(aes(x = microhabitat, y = rel_std_ovi)) +
  geom_col(aes(fill = species), position = "dodge2", color = "black") +
  scale_fill_manual(labels = c("P. napi", "P. rapae"),
                    name = "Butterfly",
                    values = c("deepskyblue",
                               "goldenrod")) +
  labs(y = "Relative<br>frequency<br>of ovipositions") +
  scale_x_discrete(labels = c("C", "SC+SO", "O"),
                   name = "Microhabitat") +
  scale_y_continuous(breaks = breaks_pretty(n = 3)) +
  theme(legend.text = element_text(face = "italic"),
        axis.title.y = element_markdown()))


# B: foliar temperature during oviposition -----------------------------------
(fol.temp <- ovi.temp %>% 
  mutate (species = recode (species,
                            PN = "P. napi",
                            PR = "P. rapae")) %>% 
  ggplot (aes (x = species, y = T_und)) +
  geom_boxplot (aes (fill = species),
                outlier.shape = NA) +
  geom_text(aes(x = 1.5, y = 40, label = "*"),
            size = 5) +
  scale_fill_manual(values = c("P. napi" = "deepskyblue",
                               "P. rapae" = "goldenrod1"),
                    name = "Butterfly") +
  labs (y = "Foliar temperature<br>(ºC)",
        x = "Butterfly species") +
  theme (text = element_text (size = 10),
         axis.text.x = element_text (size = 10,
                                     face = "italic"),
         axis.title.y = element_markdown (size = 10),
         legend.text = element_text(face = "italic"),
         panel.border = element_rect (fill = NA)) +
  scale_y_continuous(expand = expand_scale (mult = 0.1)))


# C: daily tmax at the microhabitat level ---------------------------------
(tmax.plot <- sensor %>%
  mutate(site_patch = paste(site, microhabitat, sep = "_")) %>% 
  group_by(site, microhabitat, sensor, winter_year, winter_jday) %>% 
  summarise(tmean = mean(TEMP, na.rm = T),
            tmax = max(TEMP, na.rm = T))%>% 
  mutate(microhabitat = factor(microhabitat,
                         levels = c("C", "SC", "SO", "O")),
         site = if_else(site == "Ld", "Lowland", "Mid-elevation"),
         site = factor(site,
                       levels = c("Mid-elevation", "Lowland"))) %>% 
  ungroup() %>% 
  boxplot(x = microhabitat, y = tmax, fill = site, res = F, outlier.shape = NA) +
  scale_fill_manual(aesthetics = c("fill", "color"),
                    values = c("forestgreen", "darkorchid"),
                    name = "Site") +
  labs(x = "Microhabitat",
       y = "Daily<br><i>T<sub>max</sub></i><br>(ºC)") +
  theme(axis.title.y = element_markdown()))


# D: foliar temperature ---------------------------------------------------
(fol.temp.plot <- data.leaves %>% 
  mutate(microhabitat = factor(microhabitat,
                               levels = c("C", "SC", "SO", "O")),
         Site = if_else(site == "CJ", "Mid-elevation", "Lowland"),
         Site = factor(Site, levels = c("Mid-elevation", "Lowland"))) %>% 
  boxplot(x = microhabitat, y = obv_temp, fill = Site, res = F,
          outlier.shape = NA) +
  scale_fill_manual(aesthetics = c("fill", "color"),
                    values = c("forestgreen", "darkorchid")) +
  guides(fill = "none", color = "none"))


# E: mean offset (microhabitat - macroclimate) ----------------------------
(offset.plot <- sensor %>%
  mutate(site_micro = paste(site, microhabitat, sep = "_")) %>% 
  group_by(site, microhabitat, sensor, winter_year, winter_jday) %>% 
  summarise(tmean = mean(TEMP, na.rm = T),
            tmax = max(TEMP, na.rm = T))%>% 
  left_join(meteo, by = c("site" = "Site", "winter_year" = "Year",
                          "winter_jday" = "jday")) %>% 
  mutate(mean_offset = tmean - TM,
         max_offset = tmax - TX,
         microhabitat = factor(microhabitat,
                         levels = c("C", "SC", "SO", "O")),
         site = if_else(site == "Ld", "Lowland", "Mid-elevation"),
         site = factor(site,
                       levels = c("Mid-elevation", "Lowland"))) %>% 
  ungroup() %>% 
  boxplot(x = microhabitat, y = mean_offset, fill = site, res = F,
          outlier.shape = NA) +
  geom_hline(aes(yintercept = 0)) +
  scale_fill_manual(aesthetics = c("fill", "color"),
                    values = c("forestgreen", "darkorchid"),
                    name = "Site") +
  labs(y = "Thermal<br>offset<br>(K)",
       x = "Microhabitat") +
  theme(axis.title.y = element_markdown()))



# F: offset at foliar level (leaf - air) ----------------------------------
(fol.offset.plot <- data.leaves %>% 
  left_join(data.plants) %>% 
  mutate(microhabitat = factor(microhabitat,
                               levels = c("C", "SC", "SO", "O")),
         thermal_offset = obv_temp - air_temp,
         Site = if_else(site == "CJ", "Mid-elevation", "Lowland"),
         Site = factor(Site, levels = c("Mid-elevation", "Lowland"))) %>% 
  boxplot(y = thermal_offset, x = microhabitat, fill = Site, res = F,
          outlier.shape = NA) +
  scale_fill_manual(aesthetics = c("fill", "color"),
                    values = c("forestgreen", "darkorchid")) +
  geom_hline(aes(yintercept = 0)) +
  guides(fill = "none", color = "none"))




# G: Daily variability on microhabitat thermal profiles -------------------
(sd.boxplot <- sensor %>%
  mutate(site_microhabitat = paste(site, microhabitat, sep = "_")) %>% 
  group_by(site, microhabitat, sensor, winter_year, winter_jday) %>% 
  mutate(microhabitat = factor(microhabitat,
                         levels = c("C", "SC", "SO", "O")),
         site = if_else(site == "Ld", "Lowland", "Mid-elevation"),
         site = factor(site,
                       levels = c("Mid-elevation", "Lowland"))) %>% 
  summarise(sd = sd(TEMP, na.rm = T)) %>% 
  ungroup() %>% 
  boxplot(x = microhabitat, y = sd, fill = site, res = F,
          outlier.shape = NA, ymax = 15) +
  scale_fill_manual(aesthetics = c("fill", "color"),
                    values = c("forestgreen", "darkorchid"),
                    name = "Site") +
  labs(y = "Daily SD<br>of temperature<br>(K)",
       x = "Microhabitat") +
  ylim(0, 15) +
  theme(axis.title.y = element_markdown()))


# H: daily heterogeneity of foliar temperatures of the same micros --------
(leafsd.plot <- data.leaves %>% 
  pivot_longer(cols = c(obv_temp, rev_temp),
               names_to = "side",
               values_to = "leaf_temp") %>% 
  mutate(microhabitat = factor(microhabitat,
                               levels = c("C", "SC", "SO", "O")),
         Site = if_else(site == "CJ", "Mid-elevation", "Lowland"),
         Site = factor(Site, levels = c("Mid-elevation", "Lowland"))) %>%
  group_by(Site, microhabitat, jday) %>% 
  summarise(sd_leaf = sd(leaf_temp, na.rm = T)) %>% 
  ungroup() %>% 
  boxplot(y = sd_leaf, x = microhabitat, fill = Site, res = F,
          outlier.shape = NA) +
  scale_fill_manual(aesthetics = c("fill", "color"),
                    values = c("forestgreen", "darkorchid")))


# Figure assemblage -------------------------------------------------------
(microclim <- ovi.plot + guides(fill = "none") +
  fol.temp + labs(y = "Oviposition<br>temperature<br>(ºC)") +
  tmax.plot + fol.temp.plot + labs(y = "Temperature<br>(ºC)",
                                   x = "Microhabitat") +
  offset.plot + fol.offset.plot + labs(y = "Thermal<br>offset<br>(K)",
                                       x = "Microhabitat") +
  sd.boxplot + leafsd.plot + labs(y = "Thermal<br>heterogeneity<br>(Daily SD, K)",
                                  x = "Microhabitat") +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(text = element_text (size = 10),
        plot.margin = margin (0, 10, 0, 5, "pt"),
        panel.border = element_rect (fill = NA),
        axis.title.y = element_markdown(angle = 90),
        axis.text.y = element_text (size = 10),
        plot.tag = element_text (size = 10,
                                 face = "bold"),
        plot.title = element_text(size = 12)))






