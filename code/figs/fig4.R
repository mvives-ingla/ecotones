############
## Fig. 3 ##
############

# Packages ----------------------------------------------------------------
library(tidyverse)
library(ggtext)
library(broom)
library(readxl)


# Data --------------------------------------------------------------------
## Field campaign data at the host plant level
data.plants <- read.csv("data/data_plants.csv")
## Field campaign data at the foliar level
data.leaves <- read.csv("data/data_leaves.csv")
## Thermal records from the ground
soilgrad <- read_excel ("data/feedbacks_united_2019.05.07.xlsx")


theme_set(theme_classic())

# b1: basal leaf - apical leaf --------------------------------------------
dif_ba_data <- data.leaves %>%
  filter(site == "AE") %>% 
  left_join(data.plants) %>% 
  group_by(site, microhabitat, jday) %>%
  mutate(nleaftype = n_distinct(leaf_type)) %>%
  filter (leaf_type %in% c("b", "m", "a"),
          nleaftype == 3,
          stem_length > 0) %>% 
  group_by(site, microhabitat, jday, leaf_type, ind) %>% 
  summarise(obv_temp = mean(obv_temp, na.rm = T)) %>%
  pivot_wider(names_from = leaf_type,
              values_from = obv_temp) %>% 
  mutate(bm = b - m,
         ba = b - a,
         ma =  m - a) %>% 
  pivot_longer(cols = bm:ma,
               names_to = "diff_type",
               values_to = "diff_obv_temp") %>% 
  filter(diff_type == "ba", !is.na(diff_obv_temp)) %>% 
  mutate(
    Microhabitat = factor(microhabitat,
                          levels = c("C", "SC", "SO", "O")),
    Microhabitat = if_else(microhabitat %in% c("SO", "SC"),
                           "SO+SC", microhabitat),
    Microhabitat = factor(Microhabitat,
                          levels = c("C", "SO+SC", "O"))
  ) %>% 
  filter(Microhabitat != "C")

dif_ba_test <- dif_ba_data %>% 
  lm(diff_obv_temp ~ Microhabitat, data = .) %>% 
  # tidy() %>% 
  # filter(term == "MicrohabitatO") %>% 
  glance() %>% 
  mutate(pval = if_else(round(p.value, 2)<= 0.05, "*", "")) %>% 
  select(pval) %>% 
  unlist()
  

(dif_ba_plot <- dif_ba_data %>%
  ggplot(aes(x = Microhabitat, y = diff_obv_temp, fill = Microhabitat)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(aes(yintercept = 0)) +
  geom_text(aes(x = 1.5, y = 6, label = dif_ba_test)) +
  scale_fill_manual(values = c("deepskyblue", "goldenrod")) +
  scale_y_continuous(breaks = c(0, 5), limits = c(-4, 6)) +
  labs(y = "Basal - Apical (K)") +
  guides(fill = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(color = "black")))

ggsave(filename = "figures/fig4_b1.svg",
       plot = dif_ba_plot,
       width = 5,
       height = 5,
       units = "cm",
       dpi = 600)



# b2: temperature change from the ground ----------------------------------

soilabs_model <- soilgrad %>%
  filter(site == "AE", rad == "R", wind == "NW")  %>% 
  lm (formula =  temp ~ sinh (0.1*(100 - height)), data = .) %>% 
  summary() %>% 
  glance() %>% 
  select(r.squared, p.value) %>% 
  mutate(rsq = round(r.squared, 2),
         pval = round(p.value, 4)) %>% 
  mutate (pval = if_else(pval == 0, "<0.0001", paste ("=", pval)),
          lab = paste0("*R*<sup>2</sup>=", rsq, "<br>*p*", pval))  


(soilabs_plot <- soilgrad %>%
  filter((site == "AE"))  %>% 
  unite (col = "rad_wind", rad, wind, sep = "+", remove = F) %>% 
  mutate (rad_wind = factor (rad_wind,
                             levels = c("R+NW", "NR+W"))) %>% 
  filter (!is.na (rad_wind)) %>% 
  ggplot (aes (y = temp, x = height)) +
  geom_point (aes (color = rad_wind), size = 0.7, alpha = 0.35) +
  geom_smooth (aes (color = rad_wind, fill = rad_wind),
               method = "lm",
               formula = y ~ sinh (0.1*(100 - x))) +
  scale_color_manual (aesthetics = c("colour", "fill"),
                      values = c ("goldenrod", "deepskyblue"))  +
    geom_richtext(x = 102,
                  y = 47,
                  label = soilabs_model$lab[1],
                  hjust = 1,
                  vjust = 1,
                  color = "goldenrod",
                  label.color = NA) +
  labs (y = "Temperature (ºC)",
        x = "Height (cm)") +
  guides(color = "none", fill = "none") +
  coord_flip() +
  theme (axis.title.x = element_blank(),
         axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
         panel.background = element_rect(color = "black")))

ggsave(filename = "figures/fig4_b2.svg",
       plot = soilabs_plot,
       width = 5,
       height = 5,
       units = "cm",
       dpi = 600)


# c2: leaf side vs microhabitat temperature -------------------------------

obvrev_temp <- data.leaves %>% 
  mutate(dif_temp = obv_temp - rev_temp,
         Microhabitat = factor(microhabitat,
                               levels = c("C", "SC", "SO", "O")),
         Microhabitat = if_else(Microhabitat %in% c("SO", "SC"),
                                "SO+SC", microhabitat),
         Microhabitat = factor(Microhabitat,
                               levels = c("C", "SO+SC", "O"))) %>% 
  filter(Microhabitat != "C", site == "AE",
         !is.na(dif_temp), !is.na(SENS_TX)) 

obvrev_fit <- obvrev_temp %>%
  filter(Microhabitat == "O") %>% 
  lm(dif_temp ~ poly(SENS_TX, 2), data = .) %>% 
  glance() %>% 
  mutate(rsq = round(r.squared, 2),
         pval = if_else(p.value < 0.0001,
                        "<0.0001",
                        as.character(round(p.value, 2))),
         fit = paste0("<i>R</i><sup>2</sup>=", rsq,
                      ",<br><i>p</i>", pval))


(obvrev_tmax <- obvrev_temp %>% 
  ggplot(aes(x = SENS_TX, y = dif_temp, color = Microhabitat)) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_smooth(se = F, method = "lm", formula = y ~ poly(x, 2)) +
  geom_richtext(aes(x = 19, y = 10, label = obvrev_fit$fit[1]),
                fill = NA, label.color = NA, hjust = 0,
                color = "goldenrod") +
  scale_color_manual(values = c("SO+SC" = "deepskyblue",
                                "O" = "goldenrod"),
                     name = "Microhabitat") +
  labs(x = "Microhab. Tmax (ºC)",
       y = "Upper - Under (K)") +
  guides(color = "none") +
  theme(panel.background = element_rect(color = "black")))


ggsave(filename = "figures/fig4_c2.svg",
       plot = obvrev_tmax,
       width = 5,
       height = 5,
       units = "cm",
       dpi = 600)


# c3: leaf - air ----------------------------------------------------------
leafair <- data.leaves %>% 
   left_join(data.plants) %>% 
   mutate(Microhabitat = factor(microhabitat,
                                levels = c("C", "SC", "SO", "O")),
          Microhabitat = if_else(Microhabitat %in% c("SO", "SC"),
                                 "SO+SC", microhabitat),
          Microhabitat = factor(Microhabitat,
                                levels = c("C", "SO+SC", "O")),
          thermal_offset = obv_temp - air_temp) %>% 
   filter(site == "AE", Microhabitat != "C")

leafair.test <- leafair %>% 
  lm(thermal_offset ~ Microhabitat, data = .) %>% 
  glance() %>% 
  mutate(ast = if_else(p.value <= 0.05, "*", ""))

(leafair.plot <- leafair %>% 
  ggplot(aes(y = thermal_offset, x = Microhabitat, fill = Microhabitat)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("deepskyblue", "goldenrod")) +
  geom_hline(aes(yintercept = 0)) +
  geom_text(x = 1.5, y = 15, label = leafair.test$ast[1]) +
  labs(y = "Leaf - Air (K)") +
  guides(fill = "none") +
    theme (axis.title.x = element_blank(),
           axis.text.x = element_blank(),
           axis.ticks.x = element_blank(),
           panel.background = element_rect(color = "black")))


ggsave(filename = "figures/fig4_c3.svg",
       plot = leafair.plot,
       width = 5,
       height = 5,
       units = "cm",
       dpi = 600)
