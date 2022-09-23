############
## Fig. 4 ##
############


# Packages ----------------------------------------------------------------
library(tidyverse)
library(ggtext)
library(patchwork)

# Functions ---------------------------------------------------------------
source("code/figs/boxplot_ecotones_v11.R")
source("code/figs/result_bivar_in_subtitle.R")


# Data --------------------------------------------------------------------
## Field campaign data at the host plant level
data.plants <- read.csv("data/data_plants.csv")
## Field campaign data at the foliar level
data.leaves <- read.csv("data/data_leaves.csv")


# A: soil humidity boxplot -----------------------------------------------------------
(soilhum.plot <- data.plants %>% 
  mutate(Microhabitat = factor(microhabitat,
                               levels = c("C", "SC", "SO", "O")),
         Site = if_else(site == "AE", "Lowland", "Mid-elevation"),
         Site = factor(Site, levels = c("Mid-elevation",
                                        "Lowland"))) %>% 
  boxplot(x = Microhabitat, y = soil_hum_org, fill = Site,
          outlier.shape = NA, res = F) +
   guides(color = "none")) 


# B+C: soil humidity trend ------------------------------------------------
(soilhum.trend <- data.plants %>%
  mutate(microhabitat = factor(microhabitat,
                               levels = c("C", "SC", "SO", "O"))) %>% 
  filter(microhabitat != "C") %>% 
  split(.$species) %>% 
  map(~ ggplot(data = .x, aes(x = jday, y = soil_hum_org, color = microhabitat)) +
        geom_point(size = 0.5) +
        geom_smooth(aes(linetype = microhabitat), se = F, span = 1) +
        scale_color_manual(values = c("deepskyblue", "deepskyblue", "goldenrod")) +
        scale_linetype_manual(values = c("dashed", "solid", "solid")) +
        scale_x_continuous(breaks = seq(60, 300, by = 60)) +
        labs(y = "Soil<br>humidity<br>(%)",
             x = "Ordinal day") +
        theme_classic() +
        theme(text = element_text(size = 10),
              plot.margin = margin(0, 10, 0, 5, "pt"),
              panel.border = element_rect(fill = NA),
              axis.title.x = element_text(size = 10),
              axis.title.y = element_markdown(angle = 90, size = 10),
              axis.text.y = element_text (size = 10))))

# D: foliar water content boxplot --------------------------------------------
(wat.cont.plot <- data.leaves %>%
  mutate(Microhabitat = factor(microhabitat,
                               levels = c("C", "SC", "SO", "O")),
         Site = if_else(site == "AE", "Lowland", "Mid-elevation"),
         Site = factor(Site, levels = c("Mid-elevation",
                                        "Lowland"))) %>%
  filter(phenology != "sen",
         leaf_state == "green") %>%
  boxplot(y = water_content, ymax = 12, x = Microhabitat,
          fill = Site, outlier.shape = NA, res = F) +
  ylim(1, 12) +
   guides(color = "none"))


# E+F: water content trend ------------------------------------------------
(wat.cont.trend <- data.leaves %>%
  mutate(microhabitat = factor(microhabitat,
                               levels = c("C", "SC", "SO", "O"))) %>% 
  unite ("micro_phase", microhabitat, cycle_phase, remove = F) %>% 
  filter(microhabitat != "C") %>% 
  split(.$species) %>% 
  map(~ ggplot(data = .x, aes(x = jday, y = water_content, color = microhabitat,
                              group = micro_phase)) +
        geom_point(size = 0.5) +
        geom_smooth(aes(linetype = microhabitat), se = F, span = 1) +
        scale_color_manual(values = c("deepskyblue", "deepskyblue", "goldenrod")) +
        scale_linetype_manual(values = c("dashed", "solid", "solid")) +
        scale_x_continuous(breaks = seq(60, 300, by = 60)) +
        labs(y = "Foliar water<br>content<br>(g g<sup>-1</sup> DW)",
             x = "Ordinal day") +
        theme_classic()))

# assemblage --------------------------------------------------------------

(Fig.hostplant <- soilhum.plot + labs(y = "Soil<br>humidity<br>(%)",
                      x = NULL) +
  wat.cont.plot + labs(y = "Foliar water<br>content<br>(g g<sup>-1</sup> DW)",
                       x = "Microhabitat") +
  soilhum.trend$`Alliaria petiolata` + labs(x = NULL) +
  wat.cont.trend$`Alliaria petiolata` + 
  soilhum.trend$`Lepidium draba` + labs(x = NULL) +
  wat.cont.trend$`Lepidium draba` + 
  plot_layout(ncol = 3, byrow = F, guides = "collect") +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")") &
  theme(text = element_text(size = 10),
        plot.margin = margin(0, 10, 0, 5, "pt"),
        panel.border = element_rect(fill = NA),
        plot.title = element_text(size = 10, face = "italic"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_markdown(angle = 90, size = 10),
        axis.text.y = element_text (size = 10),
        plot.tag = element_text (size = 10,
                                 face = "bold"),
        legend.position = "bottom"))

