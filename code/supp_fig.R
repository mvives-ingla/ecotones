###########################
## Supplementary figures ##
###########################


# packages ----------------------------------------------------------------
library(tidyverse)
library(ggtext)
library(colorspace)
library(patchwork)
library(scales)
library(moments)
library(zoo)
library(readxl)
library(emmeans)
library(multcomp)
library(broom)


# functions ---------------------------------------------------------------
source("code/boxplot_ecotones_v11.R")
source("code/result_bivar_in_subtitle.R")


# data --------------------------------------------------------------------
## field data at the microhabitat level
data.microh <- read.csv("data/data_microhabitat.csv")
## field data at the plant level
data.plants <- read.csv("data/data_plants.csv")
## field data at the foliar level
data.leaves <- read.csv("data/data_leaves.csv")
## reproductive output of the monitored host plants
data.fruits <- read.csv ("data/data_fruits.csv")

## thermal records in the field at different heights
soilgrad <- read_excel ("data/feedbacks_united_2019.05.07.xlsx")

## butterfly abundance in the field
but.pheno <- read.csv("data/but_pheno.csv")

## sensor data
sensor <- read.csv("data/sensors_hourly.csv") 

## correspondence of each sensor to the microhabitat type
sensor_patch <- read_delim("data/table_yearsens.txt",
                           delim = " ",
                           col_names = T) %>%
  rename(microhabitat = patch2)

## Data of simulated thermal mortality
fullmort <- read_csv("data/sim_mortalities/mort_min_tdtsim_sensors_full2021-06-25.csv") %>% 
  bind_rows(read_csv("data/sim_mortalities/mort_min_tdtsim_sensors_full2021-06-24.csv")) %>% 
  bind_rows(read_csv("data/sim_mortalities/mort_min_tdtsim_sensors_full2021-06-23.csv"))

## Simulated data
tdtdata <- read_csv("data/sim_mortalities/tdtsim_full2021-06-25.csv") %>% 
  bind_rows(read_csv("data/sim_mortalities/tdtsim_full2021-06-24.csv")) %>% 
  bind_rows(read_csv("data/sim_mortalities/tdtsim_full2021-06-23.csv"))

## Mean annual series of microclimatic field records at min resolution
sensors_min_spline <- readRDS("data/sim_mortalities/sensors_minute_spline_2021-06-25.RDS")

dayjday <- fullmort %>% 
  distinct(sensor, winter_jday) %>% 
  group_by(sensor) %>% 
  mutate(day = min_rank(winter_jday))

sensors_spline <- sensors_min_spline %>% 
  mutate(day = floor((hour-1)/24) + 1) %>% 
  left_join(dayjday)





theme_set(theme_classic())

# Fig. S3: microhabitat categories depending on vegetal cover -------------
cn.gr.data <- data.microh %>% 
  mutate (canopy = case_when ((microhabitat == "C" | microhabitat == "SC") ~ ori_mean_all,
                              TRUE ~ NA_real_),
          herb = case_when ((microhabitat == "SO" | microhabitat == "O") ~ hei_ground_cover,
                            TRUE ~ NA_real_),
          microhabitat = factor (microhabitat, levels = c ("C", "SC", "SO", "O"))) %>% 
  select (species, microhabitat, jday, herb, canopy) %>% 
  split (.$species)

(cn.gr.plot <- data.microh %>% 
  mutate (canopy = case_when ((microhabitat == "C" | microhabitat == "SC") ~ ori_mean_all,
                              TRUE ~ NA_real_),
          herb = case_when ((microhabitat == "SO" | microhabitat == "O") ~ hei_ground_cover,
                            TRUE ~ NA_real_),
          microhabitat = factor (microhabitat, levels = c ("C", "SC", "SO", "O"))) %>% 
  dplyr::select (species, microhabitat, jday, herb, canopy) %>% 
  pivot_longer(cols = herb:canopy,
               names_to = "var",
               values_to = "value") %>% 
  split(.$species) %>% 
  map(boxplot,
      x = microhabitat,
      y = value,
      res = F,
      outlier.shape = NA) %>% 
  map(~ .x +
        scale_y_continuous (limits = c(0, 100),
                            name = "Canopy cover [%]",
                            sec.axis = sec_axis (~.*1,
                                                 name = "Herb height [cm]")) +
        scale_x_discrete(name = "Microhabitat") +
        geom_vline (aes (xintercept = 2.5),
                    linetype = "dashed")))


(cn.gr.trend <- cn.gr.data %>% 
  map (~ ggplot (data = ., aes (x = jday)) +
         geom_point (aes (y = canopy, color = microhabitat),
                     size = 0.5) +
         geom_smooth (aes (y = canopy, color = microhabitat),
                      se = F) +
         geom_point (aes (y = herb, color = microhabitat),
                     size = 0.5) +
         geom_smooth (aes (y = herb, color = microhabitat),
                      se = F,
                      linetype = "dashed") +
         scale_y_continuous (limits = c(0, 100),
                             name = "Canopy cover [%]",
                             sec.axis = sec_axis (~.*1,
                                                  name = "Herb height [cm]")) +
         scale_x_continuous (name = "Julian day",
                             breaks = seq (from = 60, to = 300, by = 60)) +
         scale_color_manual (name = "Microhabitat",
                             values = heat_hcl (4,
                                                h = c(270, 30),
                                                l = 60,
                                                c = c (50, 50))) +
         labs (title = .$species) +
         theme_classic() +
         theme (axis.line.y.right = element_line (linetype = "dashed"))))


(FigS3 <- cn.gr.plot [["Alliaria petiolata"]] + cn.gr.trend [["Alliaria petiolata"]] +
  cn.gr.plot [["Lepidium draba"]] + cn.gr.trend [["Lepidium draba"]] +
  plot_layout (ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme (text = element_text (size = 10),
         plot.tag = element_text (size = 10,
                                  face = "bold"),
         plot.title = element_text (size = 10,
                                    face = "italic"),
         axis.text = element_text (size = 10),
         axis.title = element_text (size = 10)))


# Fig. S4: ovipositions during different times of the day -----------------
source("code/ovip_models.R")



# Fig. S5: trends of thermal heterogeneity in the open microhabita --------
(soiltemp <- data.microh %>% 
   left_join(data.plants) %>% 
   filter(microhabitat == "O", site == "AE", cycle_phase != "absence") %>% 
   ggplot(aes(x = jday, y = soil_surface_R)) +
   geom_point(aes(shape = cycle_phase), alpha = 0.5) +
   scale_shape_manual(values = c(1, 19),
                      labels = c("spring plants", "summer resprout"),
                      name = "Plant stage") +
   geom_smooth(se = F, color = "black") +
   labs(x = "Julian day",
        y = "Soil\ntemperature\n(ºC)"))

(microhtemp <- data.leaves %>% 
   left_join(data.plants) %>% 
   mutate(Microhabitat = factor(microhabitat,
                                levels = c("C", "SC", "SO", "O")),
          therm_off = obv_temp-air_temp,
          obv_rev = obv_temp-rev_temp) %>% 
   filter(Microhabitat == "O", species == "Lepidium draba",
          cycle_phase != "absence") %>% 
   ggplot(aes(x = jday, y = SENS_TX, shape = cycle_phase)) +
   geom_point(alpha = 0.5) +
   scale_shape_manual(values = c(1, 19),
                      labels = c("spring plants", "summer resprout"),
                      name = "Plant stage") +
   geom_smooth(se = F, aes(shape = NULL), color = "black") +
   labs(x = "Julian day",
        y = "Microhabitat<br>air <i>T<sub>max</sub></i><br>(ºC)") +
   theme(axis.title.y = element_markdown()))

## foliar temperature
th <- data.frame(th = seq(from = 35, to = 45, by = 2.5))

(leaftemp <- data.leaves %>%
    filter(microhabitat == "O",
           site == "AE",
           cycle_phase != "absence") %>%
    pivot_longer(cols = c(obv_temp, rev_temp),
                 names_to = "side",
                 values_to = "leaf_temp") %>%
    ggplot(aes(x = jday, y = leaf_temp, shape = cycle_phase)) +
    geom_hline(data = th, aes(yintercept = th, color = th)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", se = F, color = "black") +
    scale_shape_manual(values = c(1, 19),
                       labels = c("spring plants", "summer resprout"),
                       name = "Plant stage") +
    scale_color_gradientn(colours = heat.colors(50, rev = T),
                          name = "TAB\nthreshold\n(ºC)") +
    labs(x = "Julian day",
         y = "Foliar\ntemperature\n(ºC)"))

## thermal heterogeneity
(sdleaf <- data.leaves %>% 
   pivot_longer(cols = c(obv_temp, rev_temp),
                names_to = "side",
                values_to = "leaf_temp") %>% 
   group_by(species, microhabitat, jday, cycle_phase) %>% 
   summarise(sd_leaf = sd(leaf_temp, na.rm = T),
             leaf_temp = mean(leaf_temp, na.rm = T)) %>% 
   mutate(Microhabitat = factor(microhabitat,
                                levels = c("C", "SC", "SO", "O"))) %>% 
   
   ungroup() %>% 
   filter(Microhabitat == "O", species == "Lepidium draba",
          cycle_phase != "absence") %>% 
   ggplot(aes(x = jday, y = sd_leaf, shape = cycle_phase)) +
   geom_point(alpha = 0.5) +
   scale_shape_manual(values = c(1, 19),
                      labels = c("spring plants", "summer resprout"),
                      name = "Plant stage") +
   geom_smooth(se = F, color = "black", method = "lm") +
   labs(x = "Julian day",
        y = "Foliar thermal\nheterogeneity\n(Daily SD, K)"))

## thermal offset (leaf-air)
(offset <- data.leaves %>% 
    left_join(data.plants) %>% 
    mutate(Microhabitat = factor(microhabitat,
                                 levels = c("C", "SC", "SO", "O")),
           therm_off = obv_temp-air_temp) %>% 
    ungroup() %>% 
    filter(Microhabitat == "O", species == "Lepidium draba",
           cycle_phase != "absence") %>% 
    ggplot(aes(x = jday, y = therm_off, shape = cycle_phase)) +
    scale_shape_manual(values = c(1, 19),
                       labels = c("spring plants", "summer resprout"),
                       name = "Plant stage") +
    geom_point(alpha = 0.5) +
    geom_smooth(se = F, method = "lm", color = "black") +
    labs(y = "Thermal offset\n(leaf - air,\nK)",
         x = "Julian day"))

(temppattern <- soiltemp + microhtemp + leaftemp + sdleaf + offset +
    plot_layout(ncol = 1, guides = "collect") +
    plot_annotation(tag_levels = "A"))

(temppattern <- temppattern & scale_x_continuous(limits = c(70, 300),
                                 breaks = seq(90, 300, by = 60)))


# Fig. S6&S7: seasonal patterns -------------------------------------------
## Particular functions
seas_func <- function (data) {
  df <- data %>%
    mutate (season_cycle_2 = case_when (phenology == "sen" ~ "senescent",
                                        phenology == "veg" ~ "rep",
                                        phenology == "rep" ~ "rep",
                                        phenology == "sum_ros" ~ "sum_ros",
                                        phenology == "aut_ros" ~ "aut_ros"),
            season_cycle_2 = factor (season_cycle_2,
                                     levels = c ("rep",
                                                 "senescent",
                                                 "sum_ros",
                                                 "aut_ros"))) %>%
    unite ("sp_stage", species, season_cycle_2, remove = F) %>%  
    mutate (sp_stage = factor (sp_stage, levels = c ("Alliaria petiolata_rep",
                                                     "Alliaria petiolata_senescent",
                                                     "Alliaria petiolata_sum_ros",
                                                     "Alliaria petiolata_aut_ros",
                                                     "Lepidium draba_rep",
                                                     "Lepidium draba_senescent",
                                                     "Lepidium draba_sum_ros",
                                                     "Lepidium draba_aut_ros")),
            microhabitat = factor (microhabitat, levels = c ("C", "SC", "SO", "O")))
  
  return (df)
}


seas_func_2 <- function (data) {
  df <- data %>%
    mutate (season_cycle_2 = case_when (month_name == "jun" ~ "senescent",
                                        (month_name == "mar" |
                                           month_name == "feb" |
                                           month_name == "apr" |
                                           month_name == "may") ~ "rep",
                                        (month_name == "jul" |
                                           month_name == "aug") ~ "sum_ros",
                                        (season_name == "autumn") ~ "aut_ros"),
            season_cycle_2 = factor (season_cycle_2,
                                     levels = c ("rep",
                                                 "senescent",
                                                 "sum_ros",
                                                 "aut_ros"))) %>%
    unite ("sp_stage", species, season_cycle_2, remove = F) %>% 
    mutate (sp_stage = factor (sp_stage, levels = c ("Alliaria petiolata_rep",
                                                     "Alliaria petiolata_senescent",
                                                     "Alliaria petiolata_sum_ros",
                                                     "Alliaria petiolata_aut_ros",
                                                     "Lepidium draba_rep",
                                                     "Lepidium draba_senescent",
                                                     "Lepidium draba_sum_ros",
                                                     "Lepidium draba_aut_ros")),
            microhabitat = factor (microhabitat, levels = c ("C", "SC", "SO", "O")))
  
  return (df)
} 

seas.data.leaves <- seas_func (data.leaves)
seas.data.plants <- seas_func (data.plants)
sp_stage <- levels (seas.data.plants$sp_stage)

## Foliar temperature
(seas.temp <- seas_func_2 (data.leaves) %>%
  split (.$sp_stage) %>% 
  map(boxplot,
      x = microhabitat,
      y = obv_temp,
      outlier.shape = NA,
      res = F) %>% 
  map2 (sp_stage,
        ~ .x +
          labs (title = .y)))


## Foliar thermal offset

(seas.amp <- seas_func_2 (data.leaves) %>% 
  mutate (fol_amp = obv_temp - METCAT_TX) %>% 
  split (.$sp_stage) %>% 
  map(boxplot,
        x = microhabitat,
        y = fol_amp,
        outlier.shape = NA,
        res = F) %>% 
  map2 (sp_stage,
        ~ .x +
          geom_hline (aes (yintercept = 0)) +
          labs (title = .y)))


## Soil humidity
(seas.hum <- seas_func_2 (data.plants) %>% 
  split (.$sp_stage) %>% 
    map(boxplot,
        x = microhabitat,
        y = soil_hum_org,
        outlier.shape = NA,
        res = F) %>% 
  map2 (sp_stage,
        ~ .x +
          labs (title = .y) +
          scale_x_discrete (drop = F)))



## Stem length
(seas.stem <- seas.data.plants %>%
  group_by (sp_stage, microhabitat) %>% 
  filter (!is.na (sp_stage)) %>%
  mutate (rank_height = row_number (stem_length)) %>% 
  top_n (5, rank_height) %>% 
  ungroup () %>% 
  split (.$sp_stage) %>%
    map(boxplot,
        x = microhabitat,
        y = stem_length,
        outlier.shape = NA,
        res = F) %>% 
  map2 (sp_stage,
        ~ .x +
          labs (title = .y) +
          scale_x_discrete (drop = F)))



## Foliar length
(seas.length <- seas.data.leaves %>%  
  filter (length != 0) %>% 
  split (.$sp_stage) %>% 
  map(boxplot,
        x = microhabitat,
        y = length,
        outlier.shape = NA,
        res = F) %>% 
  map2 (sp_stage,
        ~ .x +
          labs (title = .y) +
          scale_x_discrete (drop = F)))


## Foliar water content
(seas.water <- seas.data.leaves %>%
  split (.$sp_stage) %>% 
    map(boxplot,
        x = microhabitat,
        y = water_content,
        outlier.shape = NA,
        res = F) %>% 
  map2 (sp_stage,
        ~ .x +
          labs (title = .y) +
          scale_x_discrete (drop = F)))


## Foliar density
(seas.dens <- seas.data.leaves %>% 
  filter (leaf_dens > 0) %>% 
  split (.$sp_stage) %>% 
  map(boxplot,
        x = microhabitat,
        y = leaf_dens,
        outlier.shape = NA,
        res = F) %>% 
  map2 (sp_stage,
        ~ .x +
          labs (title = .y) +
          scale_x_discrete (drop = F)))


## Chl
(seas.chl <- seas.data.leaves %>% 
  split (.$sp_stage) %>% 
    map(boxplot,
        x = microhabitat,
        y = chl,
        outlier.shape = NA,
        res = F) %>% 
  map2 (sp_stage,
        ~ .x +
          labs (title = .y) +
          scale_x_discrete (drop = F)))

## Fig. S6
plots.ld <- list (seas.temp %>%
                    map (~ . +
                           labs (y = "Foliar<br>temperature<br>[ºC]") +
                           scale_y_continuous (limits = c(16, 44),
                                               breaks = pretty_breaks (n = 3))),
                  seas.amp %>%
                    map (~ . +
                           labs (y = "Foliar thermal<br>amplification<br>[ºC]") +
                           scale_y_continuous (limits = c (-6, 18),
                                               breaks = pretty_breaks (n = 3))),
                  seas.hum %>%
                    map (~ . +
                           labs (y = "Soil<br>humidity<br>[%]") +
                           scale_y_continuous (limits = c (0, 67),
                                               breaks = pretty_breaks (n = 3))),
                  seas.stem %>%
                    map (~ . +
                           labs (y = "Stem<br>length<br>[cm]") +
                           scale_y_continuous (limits = c (0, 140),
                                               breaks = pretty_breaks (n = 3))),
                  seas.length %>%
                    map (~ . +
                           labs (y = "Foliar<br>length<br>[cm]") +
                           scale_y_continuous (limits = c (0, 17),
                                               breaks = pretty_breaks (n = 3))),
                  seas.water %>%
                    map (~ . +
                           labs (y = "Foliar water<br>content<br>[g g<sup>-1</sup> DW]") +
                           scale_y_continuous (limits = c (0, 12),
                                               breaks = pretty_breaks (n = 3))),
                  seas.dens %>%
                    map (~ . +
                           labs (y = "Foliar<br>density<br>[mg cm<sup>-1</sup>]") +
                           scale_y_continuous (limits = c (0, 14.7),
                                               breaks = pretty_breaks (n = 3))),
                  seas.chl %>%
                    map (~ . +
                           labs (y = "Chlorophyll<br>content<br>[SPAD]") +
                           scale_y_continuous (limits = c (0, 70),
                                               breaks = pretty_breaks (n = 3)))
)

panels.ld <- list ()

for (i in 1:length(plots.ld)) {
  n <- 4*(i-1)+1
  panels.ld[n] <- plots.ld[[i]]["Lepidium draba_rep"]
  
  n <- 4*(i-1)+2
  panels.ld[n] <- plots.ld[[i]]["Lepidium draba_senescent"]
  
  n <- 4*(i-1)+3
  panels.ld[n] <- plots.ld[[i]]["Lepidium draba_sum_ros"]
  
  n <- 4*(i-1)+4
  panels.ld[n] <- plots.ld[[i]]["Lepidium draba_aut_ros"]
  
}


panels.ld <- panels.ld %>% 
  map (~ . +
         theme (text = element_text (size = 10),
                plot.margin = margin (0, 10, 0, 5, "pt"),
                panel.border = element_rect (fill = NA),
                plot.title = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.x = element_blank (),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text (size = 10),
                plot.tag = element_text (size = 10,
                                         face = "bold"))) %>% 
  map_at (seq (from = 1, length.out = length (plots.ld), by = 4),
          ~ . +
            theme (axis.title.y = element_markdown (size = 10,
                                                    angle = 90))) %>%
  map_at ((length (panels.ld):1)[1:4],
          ~ . +
            theme (axis.title.x = element_text (size = 10),
                   axis.text.x = element_text (size = 10),
                   axis.ticks.x = element_line ()) +
            labs (x = "Microhabitat")) %>% 
  map_at (1:4,
          ~ . +
            theme (plot.title = element_text (size = 10)))

FigS6 <- panels.ld[[1]]
count <- 1
while (count < length (panels.ld)) {
  count <- count + 1
  FigS6 <- FigS6 + panels.ld[[count]]
}


(FigS6 <- FigS6 + plot_layout (ncol = 4,
                              byrow = T))


## Fig. S7: AP

plots.ap <- list (seas.temp %>%
                    map (~ . +
                           labs (y = "Foliar<br>temperature<br>[ºC]") +
                           scale_y_continuous (limits = c(13, 38),
                                               breaks = pretty_breaks (n = 3))),
                  seas.amp %>%
                    map (~ . +
                           labs (y = "Foliar thermal<br>amplification<br>[ºC]") +
                           scale_y_continuous (limits = c (-11, 8),
                                               breaks = pretty_breaks (n = 3))),
                  seas.hum %>%
                    map (~ . +
                           labs (y = "Soil<br>humidity<br>[%]") +
                           scale_y_continuous (limits = c (4, 64),
                                               breaks = pretty_breaks (n = 3))),
                  seas.stem %>%
                    map (~ . +
                           labs (y = "Stem<br>length<br>[cm]") +
                           scale_y_continuous (limits = c (0, 127),
                                               breaks = pretty_breaks (n = 3))),
                  seas.length %>%
                    map (~ . +
                           labs (y = "Foliar<br>length<br>[cm]") +
                           scale_y_continuous (limits = c (1, 22),
                                               breaks = pretty_breaks (n = 3))),
                  seas.water %>%
                    map (~ . +
                           labs (y = "Foliar water<br>content<br>[g g<sup>-1</sup> DW]") +
                           scale_y_continuous (limits = c (0, 12),
                                               breaks = pretty_breaks (n = 3))),
                  seas.dens %>%
                    map (~ . +
                           labs (y = "Foliar<br>density<br>[mg cm<sup>-1</sup>]") +
                           scale_y_continuous (limits = c (0, 25),
                                               breaks = pretty_breaks (n = 3))),
                  seas.chl %>%
                    map (~ . +
                           labs (y = "Chlorophyll<br>content<br>[SPAD]") +
                           scale_y_continuous (limits = c (0, 44),
                                               breaks = pretty_breaks (n = 3)))
)

panels.ap <- list ()
for (i in 1:length(plots.ap)) {
  n <- 4*(i-1)+1
  panels.ap[n] <- plots.ap[[i]]["Alliaria petiolata_rep"]
  
  n <- 4*(i-1)+2
  panels.ap[n] <- plots.ap[[i]]["Alliaria petiolata_senescent"]
  
  n <- 4*(i-1)+3
  panels.ap[n] <- plots.ap[[i]]["Alliaria petiolata_sum_ros"]
  
  n <- 4*(i-1)+4
  panels.ap[n] <- plots.ap[[i]]["Alliaria petiolata_aut_ros"]
  
}

panels.ap <- panels.ap %>% 
  map (~ . +
         theme (text = element_text (size = 10),
                plot.margin = margin (0, 10, 0, 5, "pt"),
                panel.border = element_rect (fill = NA),
                plot.title = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.x = element_blank (),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text (size = 10),
                plot.tag = element_text (size = 10,
                                         face = "bold"))) %>% 
  map_at (seq (from = 1, length.out = length (plots.ap), by = 4),
          ~ . +
            theme (axis.title.y = element_markdown (size = 10,
                                                    angle = 90))) %>%
  map_at ((length (panels.ap):1)[1:4],
          ~ . +
            theme (axis.title.x = element_text (size = 10),
                   axis.text.x = element_text (size = 10),
                   axis.ticks.x = element_line ()) +
            labs (x = "Microhabitat")) %>% 
  map_at (1:4,
          ~ . +
            theme (plot.title = element_text (size = 10)))

FigS7 <- panels.ap[[1]]
count <- 1
while (count < length (panels.ap)) {
  count <- count + 1
  FigS7 <- FigS7 + panels.ap[[count]]
}


(FigS7 <- FigS7 + plot_layout (ncol = 4, byrow = T) + plot_annotation (tag_levels = "A"))


FigS7



# Fig. S8: distribution, skewness and mortality --------
## A: density
sens_skew_h <- sensor %>% 
  left_join(sensor_patch) %>%
  mutate(winter_jday = if_else(winter_year == 2016, winter_jday - 1, as.double(winter_jday))) %>% #Pq el 2016 tingui la mateixa correspondència de dies que els altres anys
  group_by(site, microhabitat, Sensor_name, winter_jday, winter_hour) %>%
  summarise(ta.h = mean(TEMP, na.rm = T)) %>% 
  filter(!is.na(microhabitat)) %>% 
  group_by(site, microhabitat, Sensor_name, winter_jday) %>% 
  summarise(skew = skewness(ta.h)) %>%
  group_by(site, microhabitat) %>% 
  summarise(skew = mean(skew, na.rm = T)) %>% 
  mutate(microhabitat = factor(microhabitat, levels = c("C", "SC", "SO", "O"))) %>% 
  arrange(site, microhabitat) %>% 
  ungroup() %>% 
  mutate(x = 50,
         y = rep(seq(1, by = -0.15, length.out = 4), 2),
         skew_text = paste0("sk<sub>", microhabitat, "</sub> = ", round(skew, 2)),
         site_comp = if_else(site == "Ld", "Lowland", "Mid-elevation"),
         site_comp = factor(site_comp,
                            levels = c("Mid-elevation", "Lowland")))

(sk.dens <- sensor %>% 
  left_join(sensor_patch) %>%
  mutate(winter_jday = if_else(winter_year == 2016, winter_jday - 1, as.double(winter_jday)),
         microhabitat = factor(microhabitat,
                               levels = c("C", "SC", "SO", "O"))) %>% #Pq el 2016 tingui la mateixa correspondència de dies que els altres anys
  group_by(site, microhabitat, Sensor_name, winter_jday, winter_hour) %>%
  summarise(ta.h = mean(TEMP, na.rm = T)) %>% 
  filter(!is.na(microhabitat)) %>% 
  ggplot(data = ., aes(x = ta.h, after_stat(ndensity))) +
  geom_density(aes(color = microhabitat), alpha = 0.5) +
  geom_richtext(data = sens_skew_h,
                aes(x = x, y = y, label = skew_text, color = microhabitat),
                hjust = 1, label.color = NA, fill = NA, vjust = 0.75,
                size = 3) +
  scale_color_discrete_sequential (palette = "Viridis",
                                   l1 = 25,
                                   l2 = 80,
                                   c1 = 40,
                                   c2 = 85,
                                   rev = F,
                                   name = "Microhabitat",
                                   guide = guide_legend(
                                     override.aes = list(
                                       label = "-",
                                       linetype = 0,
                                       size = 5
                                     ))) +
  facet_wrap(vars(site_comp), ncol = 2) +
  labs(y = "Scaled density",
       x = "Temperature (ºC)"))

## B: regression of skewness (calculated from data at hourly resolution) vs mean thermal mortality
sens_sk <- sensor %>% 
  left_join(sensor_patch) %>%
  mutate(winter_jday = if_else(winter_year == 2016, winter_jday - 1, as.double(winter_jday))) %>% #Pq el 2016 tingui la mateixa correspondència de dies que els altres anys
  group_by(site, microhabitat, Sensor_name, winter_jday, winter_hour) %>%
  summarise(ta.h = mean(TEMP, na.rm = T)) %>% 
  filter(!is.na(microhabitat)) %>% 
  group_by(site, microhabitat, Sensor_name, winter_jday) %>% 
  summarise(skew = skewness(ta.h)) %>%
  group_by(site, microhabitat, Sensor_name) %>% 
  summarise(skew = mean(skew, na.rm = T)) 


sk.reg.fit <- fullmort %>% 
  group_by(sensor, sp_comp, thres, simulation) %>% 
  arrange(winter_jday, .by_group = T) %>% 
  mutate(prob_surv = surv/100,
         log_prob = log10(prob_surv),
         log_prob = if_else(log_prob == -Inf, -999, log_prob),
         log_prob_30 = rollsum(log_prob, k = 30, fill = NA, align = "left"),
         prob_30 = 10^log_prob_30,
         mort_30 = 1-prob_30,
         high_mort_5 = if_else(mort_30 >= 0.05, 1, 0),
         high_mort_1 = if_else(mort_30 >= 0.01, 1, 0),
         mort_30_filt = if_else(high_mort_1 == 1, mort_30, NA_real_)) %>% 
  filter(!is.na(prob_30)) %>% 
  group_by(sensor, sp_comp, thres, simulation) %>% 
  summarise(max_mort_30 = max(mort_30),
            tmax = max(tmax),
            days_highmort_1 = sum(high_mort_1),
            mean_mort_1 = mean(mort_30_filt, na.rm = T)) %>% 
  filter(thres == 100) %>% 
  pivot_wider(id_cols = c(sensor, thres, simulation),
              names_from = sp_comp,
              values_from = mean_mort_1)  %>% 
  rename(PN = "P. napi", PR = "P. rapae") %>% 
  mutate(dif_mort = PN - PR,
         rat_mort = PN/PR) %>% 
  left_join(sensor_patch, by = c("sensor" = "Sensor_name")) %>% 
  left_join(sens_sk, by = c("sensor" = "Sensor_name")) %>% 
  lm(rat_mort ~ skew, data = .) %>% 
  summary()

(sk.reg <- fullmort %>% 
  group_by(sensor, sp_comp, thres, simulation) %>% 
  arrange(winter_jday, .by_group = T) %>% 
  mutate(prob_surv = surv/100,
         log_prob = log10(prob_surv),
         log_prob = if_else(log_prob == -Inf, -999, log_prob),
         log_prob_30 = rollsum(log_prob, k = 30, fill = NA, align = "left"),
         prob_30 = 10^log_prob_30,
         mort_30 = 1-prob_30,
         high_mort_5 = if_else(mort_30 >= 0.05, 1, 0),
         high_mort_1 = if_else(mort_30 >= 0.01, 1, 0),
         mort_30_filt = if_else(high_mort_1 == 1, mort_30, NA_real_)) %>% 
  filter(!is.na(prob_30)) %>% 
  group_by(sensor, sp_comp, thres, simulation) %>% 
  summarise(max_mort_30 = max(mort_30),
            tmax = max(tmax),
            days_highmort_1 = sum(high_mort_1),
            mean_mort_1 = mean(mort_30_filt, na.rm = T)) %>% 
  filter(thres == 100) %>% 
  pivot_wider(id_cols = c(sensor, thres, simulation),
              names_from = sp_comp,
              values_from = mean_mort_1)  %>% 
  rename(PN = "P. napi", PR = "P. rapae") %>% 
  mutate(dif_mort = PN - PR,
         rat_mort = PN/PR) %>% 
  left_join(sensor_patch, by = c("sensor" = "Sensor_name")) %>% 
  left_join(sens_sk, by = c("sensor" = "Sensor_name")) %>% 
  ggplot(aes(x = skew, y = rat_mort)) +
  geom_point(size = 0.25, alpha = 0.5,
             position = position_jitter()) +
  geom_smooth(method = "lm") +
  geom_hline(aes(yintercept = 1)) +
  labs(y = "Ratio of mean\nthermal mortality\n(Pn/Pr)",
       x = "Skewness of microhabitat temperature"))

## assemblage
(skewsd.plot <- (sk.dens + guides(color = "none")) /
  sk.reg +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(plot.title = element_text(size = 11),
        plot.subtitle = element_markdown(size = 11),
        strip.text = element_markdown(size = 11),
        axis.title = element_text(size = 11)))



# Fig. S9: details of thermal heterogeneity -------------------------------
pal <- sequential_hcl(n = 2,
                      palette = "Viridis",
                      l1 = 25, l2 = 80, c1 = 40, c2 = 85, 
                      rev = F)

## A, B: boxplots of temperatures at different heights
plant.mean <- data.plants %>% 
  group_by(site, microhabitat, jday) %>% 
  summarise(across(.cols = c(soil_temp_couple, air_temp),
                   .fns = mean, na.rm = T))
leaf.mean <- data.leaves %>% 
  group_by(site, microhabitat, jday) %>% 
  summarise(obv_temp = mean(obv_temp, na.rm = T))

boxdata.678 <- data.microh %>% 
  left_join(leaf.mean) %>% 
  left_join(plant.mean) %>% 
  pivot_longer(cols = c(soil_temp_couple,
                        soil_surface_R,
                        soil_surface_NR,
                        obv_temp,
                        SENS_TM2H,
                        air_temp,
                        METCAT_TX),
               names_to = "measure",
               values_to = "temp") %>% 
  filter(site == "AE",
         month >= 6, month <= 8,
         ((microhabitat == "C" & measure != "soil_surface_R") |
         (microhabitat == "O" & measure != "soil_surface_NR"))) %>% 
  mutate (measure = factor (measure, levels = c ("soil_temp_couple",
                                                 "soil_surface_NR",
                                                 "soil_surface_R",
                                                 "obv_temp",
                                                 "SENS_TM2H",
                                                 "air_temp",
                                                 "METCAT_TX"))) %>%
  split (.$microhabitat)


(boxplot.678 <- boxdata.678 %>% 
  map(boxplot, x = measure, y = temp, res = F,
      outlier.shape = NA, ymax = 55)) %>% 
  map2 (pal,
        ~ .x + 
          geom_boxplot (outlier.shape = NA, fill = .y) +
         geom_vline (aes (xintercept = 1.5)) +
         labs (x = NULL,
               y = "Temperature (ºC)",
               title = NULL) +
          ylim(20, 55) +
         scale_x_discrete (labels = c (soil_temp_couple = "Soil [-10 cm]",
                                       soil_surface_NR = "Soil surface",
                                       soil_surface_R = "Soil surface",
                                       SENS_TM2H = "Air [25 cm]",
                                       obv_temp = "Leaf upperside",
                                       air_temp = "Air [100 cm]",
                                       METCAT_TX = "Air [150 cm]")) +
         theme(axis.text.x = element_text(angle = 30, hjust = 1)))

## C: thermal gradient from the ground

soilabs.model <- soilgrad %>%
  filter(site == "AE")  %>% 
  unite (col = "rad_wind", rad, wind, sep = "+", remove = F) %>% 
  mutate (rad_wind = factor (rad_wind,
                             levels = c("R+NW", "NR+W"))) %>% 
  filter (!is.na (rad_wind)) %>%
  split (.$rad_wind) %>% 
  map (~ lm (formula =  temp ~ sinh (0.1*(100 - height)), data = .)) %>% 
  map(summary) %>% 
  map(~ keep(., names(.) %in% c("fstatistic", "r.squared"))) %>%
  map(flatten_dfc) %>%
  map(~ transmute(.,
                  rsq = r.squared,
                  pval = pf(q = value, df1 = numdf, df2 = dendf,
                            lower.tail = F))) %>%
  map(~ map2_dfc (.,
                  c(2, 4),
                  ~ round(.x, .y))) %>%
  bind_rows(.id = "group") %>% 
  mutate (pval = if_else(pval == 0, "< 0.0001", paste ("=", pval)),
          lab = paste0(group, ": *R*<sup>2</sup> = ", rsq, ", *p* ", pval))  


fit.compl <- str_c(soilabs.model$lab, collapse = "<br>")

(soilabs.plot <- soilgrad %>%
  filter((site == "AE"))  %>% 
  unite (col = "rad_wind", rad, wind, sep = "+", remove = F) %>% 
  mutate (rad_wind = factor (rad_wind,
                             levels = c("R+NW", "NR+W"))) %>% 
  filter (!is.na (rad_wind)) %>% 
  ggplot (aes (x = height, y = temp)) +
  geom_point (aes (color = rad_wind), size = 0.7, alpha = 0.35) +
  geom_smooth (aes (color = rad_wind, fill = rad_wind),
               method = "lm",
               formula = y ~ sinh (0.1*(100 - x))) +
  scale_color_manual (aesthetics = c("colour", "fill"),
                      values = c (pal[2], "darkorange3"),
                      guide = guide_legend(override.aes = list(fill = NA)))  +
  labs (y = "Temperature (ºC)",
        x = "Height from the ground (cm)",
        subtitle = fit.compl) +
  theme (legend.title = element_blank (),
         legend.position = "bottom",
         plot.subtitle = element_markdown(size = 10)))

## D: differences between apical and basal leaves
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
  filter(diff_type == "ba", !is.na(diff_obv_temp))  %>% 
  mutate(Microhabitat = factor(microhabitat,
                               levels = c("C", "SC", "SO", "O")))

dif_ba_tukey <- dif_ba_data %>% 
  filter(site == "AE") %>% 
  lm(diff_obv_temp ~ Microhabitat, data = .) %>% 
  emmeans(~ Microhabitat) %>% 
  multcomp::cld(method = "Tukey", Letters = letters)


(dif_ba_plot <- dif_ba_data %>%
  ggplot(aes(x = Microhabitat, y = diff_obv_temp)) +
  geom_boxplot(aes(fill = Microhabitat), outlier.alpha = 0.5, outlier.size = 0.5) +
  geom_hline(aes(yintercept = 0)) +
  geom_text(data = dif_ba_tukey, aes(x = Microhabitat,
                                     y = 10,
                                     label = .group)) +
  scale_fill_discrete_sequential (palette = "Viridis",
                                  l1 = 25,
                                  l2 = 80,
                                  c1 = 40,
                                  c2 = 85,
                                  rev = F,
                                  name = "Microhabitat") +
  labs(y = "T<sub>basal leaf</sub> - T<sub>apical leaf</sub> (K)") +
  theme(strip.text = element_text(face = "italic"),
        axis.title.y = element_markdown()))

## E: differences between upperside and under side
obvrev_tukey <- data.leaves %>% 
  mutate(dif_temp = obv_temp - rev_temp) %>% 
  filter(species == "Lepidium draba") %>% 
  mutate(Microhabitat = factor(microhabitat,
                               levels = c("C", "SC", "SO", "O"))) %>% 
  lm(dif_temp ~ Microhabitat, data = .) %>% 
  emmeans(~ Microhabitat) %>% 
  multcomp::cld(method = "tukey", Letters = letters)


(obvrev_boxplot <- data.leaves %>% 
  mutate(dif_temp = obv_temp - rev_temp) %>% 
  filter(species == "Lepidium draba") %>% 
  mutate(
    Microhabitat = factor(microhabitat,
                          levels = c("C", "SC", "SO", "O"))) %>% 
  ggplot(aes(x = Microhabitat, y = dif_temp)) +
  geom_boxplot(aes(fill = Microhabitat),
               outlier.size = 0.5, outlier.alpha = 0.5) +
  scale_fill_discrete_sequential (palette = "Viridis",
                                  l1 = 25,
                                  l2 = 80,
                                  c1 = 40,
                                  c2 = 85,
                                  rev = F,
                                  name = "Microhabitat") +
  geom_text(data = obvrev_tukey, aes(x = Microhabitat, y = 10, label = .group)) +
  geom_hline(aes(yintercept = 0)) +
  labs(y = "T<sub>upper side</sub> - T<sub>underside</sub> (K)") +
  theme_classic() +
  theme(axis.title.y = element_markdown()))

## F: thermal difference between upper and unders vs tmax
obvrev_tmax_fit <- data.leaves %>% 
  mutate(dif_temp = obv_temp - rev_temp) %>% 
  filter(species == "Lepidium draba") %>% 
  mutate(Microhabitat = factor(microhabitat,
                               levels = c("C", "SC", "SO", "O"))) %>% 
  filter(!is.na(dif_temp),
         !is.na(SENS_TX)) %>%
  nest(data = -Microhabitat) %>% 
  mutate(lm = map(data, ~lm(dif_temp ~ poly(SENS_TX, 2), data = .)),
         fit = map(lm, glance),
         fit = map(fit, ~mutate(.,
                                rsq = round(r.squared, 2),
                                pval = if_else(p.value < 0.0001,
                                               "< 0.0001",
                                               as.character(round(p.value, 2))),
                                fit = paste0("<i>R</i><sup>2</sup> = ", rsq,
                                             ", <i>p</i> ", pval))),
         fit = map(fit, dplyr::select, fit),
         fit = map(fit, as.character)) %>% 
  unnest(fit) %>% 
  dplyr::select(Microhabitat, fit) %>% 
  mutate(fit = paste0(Microhabitat, ": ", fit),
         x = 15,
         y = c(-2.5, 5, 7.5, 10)) %>% 
  filter(Microhabitat == "O")


(obvrev_tmax <- data.leaves %>% 
  mutate(dif_temp = obv_temp - rev_temp) %>% 
  filter(species == "Lepidium draba") %>% 
  mutate(Microhabitat = factor(microhabitat,
                               levels = c("C", "SC", "SO", "O"))) %>% 
  filter(!is.na(dif_temp),
         !is.na(SENS_TX)) %>% 
  ggplot(aes(x = SENS_TX, y = dif_temp, color = Microhabitat)) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_smooth(se = F, method = "lm", formula = y ~ poly(x, 2)) +
  scale_color_discrete_sequential (palette = "Viridis",
                                   l1 = 25,
                                   l2 = 80,
                                   c1 = 40,
                                   c2 = 85,
                                   rev = F,
                                   name = "Microhabitat") +
  labs(x = "Microhabitat T<sub>max</sub> (ºC)",
       y = "T<sub>upper side</sub> - T<sub>underside</sub> (K)",
       subtitle = obvrev_tmax_fit$fit) +
  theme(axis.title.y = element_markdown(),
        axis.title.x = element_markdown(),
        plot.subtitle = element_markdown()))



## assemblage
(leaf_mosaic <- boxplot.678[[1]] + labs(title = "Closed microhabitat") +
  boxplot.678[[2]] + labs(title = "Open microhabitat") +
  soilabs.plot + labs(title = "Open microhabitat") +
  dif_ba_plot + guides(fill = "none") + theme(axis.title.y =
                                                element_markdown()) +
  obvrev_boxplot + guides(fill = "none") +
  obvrev_tmax + 
  plot_annotation(tag_levels = "A") +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom",
        plot.title = element_text(size = 10)))


# Fig. S10: trends in the host plant variables ---------------------------------
(soilsurf.trend <- data.microh %>% 
  mutate(microhabitat = factor(microhabitat,
                               levels = c("C", "SC", "SO", "O"))) %>% 
  filter(microhabitat != "C") %>% 
  split(.$species) %>% 
  map(~ ggplot(data = .x, aes(x = jday, y = soil_surface, color = microhabitat)) +
        # geom_point(size = 0.5) +
        geom_smooth(aes(linetype = microhabitat), se = F, span = 0.75) +
        scale_color_manual(values = c("deepskyblue", "deepskyblue", "goldenrod")) +
        scale_linetype_manual(values = c("dashed", "solid", "solid")) +
        scale_x_continuous(breaks = seq(60, 300, by = 60),
                           limits = c(60, 300)) +
        labs(y = "Soil surface<br>temperature<br>(ºC)",
             x = "Julian day") +
        theme(text = element_text(size = 10),
              plot.margin = margin(0, 10, 0, 5, "pt"),
              panel.border = element_rect(fill = NA),
              axis.title.x = element_text(size = 10),
              axis.title.y = element_markdown(angle = 90, size = 10),
              axis.text.y = element_text (size = 10))))

(leaftemp.trend <- data.leaves %>% 
  mutate(microhabitat = factor(microhabitat,
                               levels = c("C", "SC", "SO", "O"))) %>% 
  filter(microhabitat != "C") %>% 
  split(.$species) %>% 
  map(~ ggplot(data = .x, aes(x = jday, y = obv_temp, color = microhabitat)) +
        # geom_point(size = 0.5) +
        geom_smooth(aes(linetype = microhabitat), se = F, span = 0.75) +
        scale_color_manual(values = c("deepskyblue", "deepskyblue", "goldenrod")) +
        scale_linetype_manual(values = c("dashed", "solid", "solid")) +
        scale_x_continuous(breaks = seq(60, 300, by = 60),
                           limits = c(60, 300)) +
        labs(y = "Foliar<br>temperature<br>(ºC)",
             x = "Julian day") +
        theme(text = element_text(size = 10),
              plot.margin = margin(0, 10, 0, 5, "pt"),
              panel.border = element_rect(fill = NA),
              axis.title.x = element_text(size = 10),
              axis.title.y = element_markdown(angle = 90, size = 10),
              axis.text.y = element_text (size = 10))))



(pheno.trend <- data.microh %>% 
  dplyr::select(jday, species, microhabitat, n_veg:n_sum_ros) %>% 
  gather(key = "pheno", value = "prop", n_veg:n_sum_ros, factor_key = T) %>% 
  group_by(jday, species, microhabitat, pheno) %>% 
  summarise(prop = mean(prop, na.rm = T)) %>% 
  unite("sp_pheno",
        species,
        pheno,
        remove = F) %>% 
  filter(pheno == "n_rep") %>% 
  mutate(microhabitat = factor(microhabitat,
                               levels = c("C", "SC","SO", "O"))) %>% 
  filter(microhabitat != "C") %>% 
  split(.$sp_pheno) %>% 
  map2(c("Alliaria petiolata", "Lepidium draba"),
       ~ ggplot(data = .x, aes(x = jday, y = prop, color = microhabitat)) +
         # geom_point(size = 0.5) +
         geom_smooth(aes(linetype = microhabitat),
                     method = "loess",
                     span = 0.65,
                     se = F) +
         scale_color_manual(values = c("deepskyblue", "deepskyblue", "goldenrod")) +
         scale_linetype_manual(values = c("dashed", "solid", "solid")) +
         scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.5)) +
         scale_x_continuous(breaks = seq(60, 300, by = 60),
                            limits = c(60, 300)) +
         labs(y = "Proportion of<br>reproductive<br>individuals",
              x = "Julian day") +
         theme(text = element_text(size = 10),
               plot.margin = margin(0, 10, 0, 5, "pt"),
               panel.border = element_rect(fill = NA),
               axis.title.x = element_text(size = 10),
               axis.title.y = element_markdown(angle = 90, size = 10),
               axis.text.y = element_text (size = 10))))



(stemlength.trend <- data.plants %>% 
  mutate(microhabitat = factor(microhabitat,
                               levels = c("C", "SC", "SO", "O"))) %>% 
  unite ("micro_phase", microhabitat, cycle_phase, remove = F) %>% 
  filter(microhabitat != "C",
         cycle_phase != "absence") %>% 
  split(.$species) %>% 
  map(~ ggplot(data = .x, aes(x = jday, y = stem_length, color = microhabitat,
                              group = micro_phase)) +
        # geom_point(size = 0.5) +
        geom_smooth(aes(linetype = microhabitat), se = F, span = 0.75) +
        scale_color_manual(values = c("deepskyblue", "deepskyblue", "goldenrod")) +
        scale_linetype_manual(values = c("dashed", "solid", "solid")) +
        scale_x_continuous(breaks = seq(60, 300, by = 60),
                           limits = c(60, 300)) +
        labs(y = "Stem<br>length<br>(cm)",
             x = "Julian day") +
        theme(text = element_text(size = 10),
              plot.margin = margin(0, 10, 0, 5, "pt"),
              panel.border = element_rect(fill = NA),
              axis.title.x = element_text(size = 10),
              axis.title.y = element_markdown(angle = 90, size = 10),
              axis.text.y = element_text (size = 10))))

(chl.trend <- data.leaves %>%
  mutate(microhabitat = factor(microhabitat,
                               levels = c("C", "SC", "SO", "O"))) %>% 
  unite ("micro_phase", microhabitat, cycle_phase, remove = F) %>% 
  filter(microhabitat != "C",
         cycle_phase != "absence") %>% 
  split(.$species) %>% 
  map(~ ggplot(data = ., aes(x = jday, y = chl, color = microhabitat,
                             group = micro_phase)) +
        # geom_point(size = 0.5) +
        geom_smooth(aes(linetype = microhabitat), se = F, span = 1) +
        scale_color_manual(values = c("deepskyblue", "deepskyblue", "goldenrod")) +
        scale_linetype_manual(values = c("dashed", "solid", "solid")) +
        scale_x_continuous(breaks = seq(60, 300, by = 60),
                           limits = c(60, 300)) +
        labs(y = "Chlorophyll<br>content<br>(SPAD)",
             x = "Julian day") +
        theme(text = element_text(size = 10),
              plot.margin = margin(0, 10, 0, 5, "pt"),
              panel.border = element_rect(fill = NA),
              axis.title.x = element_text(size = 10),
              axis.title.y = element_markdown(angle = 90, size = 10),
              axis.text.y = element_text (size = 10))))



(but.plot <- but.pheno %>% 
  split (.$site) %>% 
  map(~ ggplot (data = .,
                aes (x = jday, y = Ab_index)) +
        geom_smooth (aes (color = IDesp),
                     span = 0.25,
                     se = F) +
        labs (x = "Julian day",
              y = "Abundance<br>index<br>(counts km<sup>-1</sup>)") +
        scale_color_manual(values = c("deepskyblue", "goldenrod"),
                           labels = c("<i>P. napi</i>", "<i>P. rapae</i>"),
                           name = "Butterfly") +
        scale_y_continuous (breaks = pretty_breaks (n = 3),
                            limits = c (0, 14),
                            expand = expand_scale (mult = 0.1)) +
        scale_x_continuous(breaks = seq(60, 300, by = 60),
                           limits = c(60, 300)) +
        theme(text = element_text(size = 10),
              plot.margin = margin(0, 10, 0, 5, "pt"),
              panel.border = element_rect(fill = NA),
              axis.title.x = element_text(size = 10),
              axis.title.y = element_markdown(angle = 90, size = 10),
              axis.text.y = element_text (size = 10),
              legend.text = element_markdown(size = 10))))


(Fig.trends <- soilsurf.trend$`Alliaria petiolata` + labs(title = "Mid-elevation site") +
  soilsurf.trend$`Lepidium draba` + labs(title = "Lowland site",
                                         y = NULL) +
  leaftemp.trend$`Alliaria petiolata`+
  leaftemp.trend$`Lepidium draba` + labs(y = NULL) +
  pheno.trend$`Alliaria petiolata_n_rep`+
  pheno.trend$`Lepidium draba_n_rep`+ labs(y = NULL) +
  stemlength.trend$`Alliaria petiolata`+
  stemlength.trend$`Lepidium draba`+ labs(y = NULL) +
  chl.trend$`Alliaria petiolata`+
  chl.trend$`Lepidium draba`+ labs(y = NULL) +
  but.plot$CJ +
  but.plot$AE + labs(y = NULL) +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(plot.title = element_text (size = 10),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()))

Fig.trends[[11]] <- Fig.trends[[11]] + theme(axis.title.x = element_text(size = 10),
                                             axis.text.x = element_text(size = 10))

Fig.trends[[12]] <- Fig.trends[[12]] + theme(axis.title.x = element_text(size = 10),
                                             axis.text.x = element_text(size = 10))



# Fig. S11: host-plant fruits ---------------------------------------------

(fruits.plot <- data.fruits %>%  
   mutate(microhabitat = factor(microhabitat,
                                levels = c("C", "SC", "SO", "O"))) %>% 
   filter(n_fruits < 500) %>% 
  split (.$species) %>% 
  map(boxplot,
      x = microhabitat,
      y = "n_fruits",
      outlier.shape = NA))

(FigS11 <- fruits.plot [["Alliaria petiolata"]] + labs (y = "Number of<br>fruits",
                                                      x = NULL,
                                                      title = "Alliaria petiolata") +
  fruits.plot [["Lepidium draba"]]  + labs (y = "Number of<br>fruits",
                                            x = NULL,
                                            title = "Lepidium draba") +
  plot_annotation (tag_levels = "A") +
  plot_layout (ncol = 2) &
  theme (text = element_text (size = 10),
         panel.border = element_rect (fill = NA),
         plot.margin = margin (5, 20, 5, 5, "pt"),
         plot.title = element_text (size = 10,
                                    face = "italic"),
         axis.title.y = element_markdown (size = 10,
                                          angle = 90),
         axis.text = element_text (size = 10),
         plot.tag = element_text (size = 10, 
                                  face = "bold")))


# Fig. S12: Phenology -----------------------------------------------------

(plant.pheno <- data.microh %>%
  dplyr::select (jday, species, microhabitat, n_veg:n_sum_ros) %>% 
  pivot_longer(cols = n_veg:n_sum_ros,
               values_to = "prop",
               names_to = "pheno") %>% 
  filter (!is.na (microhabitat)) %>%
  group_by (jday, species, microhabitat, pheno) %>% 
  summarise (prop = mean (prop, na.rm = T)) %>%
  ungroup () %>% 
  unite ("sp_patch",
         species,
         microhabitat,
         remove = F,
         sep = " - ") %>% 
  mutate (sp_patch = factor (sp_patch,
                             levels = c ("Alliaria petiolata - C",
                                         "Alliaria petiolata - SC",
                                         "Alliaria petiolata - SO",
                                         "Alliaria petiolata - O",
                                         "Lepidium draba - C",
                                         "Lepidium draba - SC",
                                         "Lepidium draba - SO",
                                         "Lepidium draba - O")),
          long_micro = case_when(microhabitat == "O" ~ "Open microhabitat",
                                 microhabitat == "SO" ~ "Semi-open microhabitat",
                                 microhabitat == "SC" ~ "Semi-closed microhabitat",
                                 T ~ "Closed microhabitat")) %>% 
  split (.$sp_patch) %>% 
  map (~ ggplot (data = ., aes (x = jday, y = prop)) +
         geom_point ( aes (color = pheno)) +
         geom_smooth (aes (color = pheno),
                      method = "loess",
                      span = 0.5,
                      se = F) +
         labs (title = .$long_micro,
               y = "Proportion of<br>individuals",
               x = "Julian day") +
         scale_y_continuous (breaks = seq (from = 0, to = 1, by = 0.5)) +
         scale_x_continuous(breaks = seq(60, 300, by = 60),
                            limits = c(60, 300)) +
         theme (text = element_text(size = 10),
                plot.margin = margin(0, 10, 0, 5, "pt"),
                plot.title = element_text(size = 10),
                panel.border = element_rect(fill = NA),
                axis.title.x = element_text(size = 10),
                axis.title.y = element_markdown(angle = 90, size = 10),
                axis.text.y = element_text (size = 10),
                legend.text = element_markdown(size = 10)) +
         scale_color_manual (values = heat_hcl (4,
                                                h = c (120, 525),
                                                l = 50,
                                                c = c (80, 80)),
                             labels = c ("Ros.", "Rep", "Sen.", "Summer\nros."),
                             name = "Phenology")))

(but.plot <- but.pheno %>% 
    split (.$site) %>% 
    map(~ ggplot (data = .,
                  aes (x = jday, y = Ab_index)) +
          geom_smooth (aes (color = IDesp),
                       span = 0.25,
                       se = F) +
          geom_point(aes(color = IDesp)) +
          labs (x = "Julian day",
                y = "Abundance<br>index<br>(counts km<sup>-1</sup>)") +
          scale_color_grey(labels = c("<i>P. napi</i>", "<i>P. rapae</i>"),
                             name = "Butterfly") +
          scale_y_continuous (breaks = pretty_breaks (n = 3),
                              limits = c (0, 14),
                              expand = expand_scale (mult = 0.1)) +
          scale_x_continuous(breaks = seq(60, 300, by = 60),
                             limits = c(60, 300)) +
          theme(text = element_text(size = 10),
                plot.margin = margin(0, 10, 0, 5, "pt"),
                panel.border = element_rect(fill = NA),
                axis.title.x = element_text(size = 10),
                axis.title.y = element_markdown(angle = 90, size = 10),
                axis.text.y = element_text (size = 10),
                legend.text = element_markdown(size = 10))))

## assemblage
(FigS12 <- plant.pheno$`Alliaria petiolata - C` +
  plant.pheno$`Lepidium draba - C` +
    theme(axis.title.y = element_blank()) +
  plant.pheno$`Alliaria petiolata - SC` +
  plant.pheno$`Lepidium draba - SC` +
    theme(axis.title.y = element_blank()) +
  plant.pheno$`Alliaria petiolata - SO` +
  plant.pheno$`Lepidium draba - SO` +
    theme(axis.title.y = element_blank()) +
  plant.pheno$`Alliaria petiolata - O` +
  plant.pheno$`Lepidium draba - O` +
    theme(axis.title.y = element_blank()) +
  but.plot$CJ +
  but.plot$AE +
    theme(axis.title.y = element_blank()) +
  plot_layout(nrow = 5,
              byrow = T,
              guides = "collect") +
  plot_annotation(tag_levels = "A"))



# Fig. S13: summary of the z and CTmax values of the simulations ----------
tdtmeans <- tdtdata %>% 
  group_by(sp) %>% 
  summarise(across(c(ctmax, z), .fns = mean)) %>% 
  mutate(Data = "sim")

tdtreal <- data.frame(sp = c("PN", "PR"),
                      z = c(4.10, 5.10),
                      ctmax = c(51.08, 53.48),
                      Data = c("real", "real")) %>% 
  add_row(tdtmeans)

(simz <- tdtdata %>% 
  ggplot(aes(x = z)) +
  geom_density(aes(fill = sp), alpha = 0.5) +
  geom_vline(data = tdtreal, aes(xintercept = z, color = sp, linetype = Data)) +
  scale_fill_manual(aesthetics = c("fill", "color"),
                    values = c("PN" = "deepskyblue",
                               "PR" = "goldenrod1"),
                    name = "Species",
                    labels = c("<i>P. napi</i>", "<i>P. rapae</i>")) +
  labs(x = "<i>z</i>",
       y = "Denisty") +
  guides(color = guide_legend(override.aes = list(linetype = 0))) +
  theme(axis.title.x = element_markdown(),
        legend.text = element_markdown()))

(simct <- tdtdata %>% 
  ggplot(aes(x = ctmax)) +
  geom_density(aes(fill = sp), alpha = 0.5) +
  geom_vline(data = tdtreal, aes(xintercept = ctmax, color = sp, linetype = Data)) +
  scale_fill_manual(aesthetics = c("fill", "color"),
                    values = c("PN" = "deepskyblue",
                               "PR" = "goldenrod1"),
                    name = "Species",
                    labels = c("<i>P. napi</i>", "<i>P. rapae</i>")) +
  labs(x = "CT<sub>max</sub>",
       y = "Denisty") +
  guides(color = guide_legend(override.aes = list(linetype = 0))) +
  theme(axis.title.x = element_markdown(),
        legend.text = element_markdown()))

(sims <- simz + simct + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A"))


# Fig. S14: examples of acute vs chronic stresses ------------------------------
(examp.temp <- sensors_spline %>% 
  mutate(daily_hour = hour - (day-1)*24 - 1) %>%
  filter(winter_jday %in% c(168, 217),
         sensor %in% c("ld_caseta_semioberta", "LD_tossal")) %>% 
  mutate(sensor = case_when(sensor == "LD_tossal" ~ "s9 - Ld - O",
                            TRUE ~ "s10 - Ld - O"),
         sensor = factor(sensor,
                         levels = c("s9 - Ld - O", "s10 - Ld - O"))) %>% 
  split(.$winter_jday) %>% 
  map2(c("17 June", "5 Aug"),
       ~ ggplot() +
         geom_tile(data = data.frame(y = seq(37.5, 45, by = 0.1)),
                   aes(x = 12,
                       y = y,
                       height = 0.1,
                       width = 24,
                       fill = y)) +
         scale_fill_gradientn(colours = heat.colors(50, rev = T),
                              name = "Temperature (ºC)",
                              guide = guide_colorbar(direction = "horizontal",
                                                     title.position = "top")) +
         geom_line(data = .x, aes(x = daily_hour, y = ta.min,
                                  linetype = sensor)) +
         geom_hline(aes(yintercept = 41.2)) +
         labs(x = "Hour of the day",
              y = "Temperature (ºC)",
              title = .y) +
         guides(linetype = guide_legend(title = "Sensor")) +
         scale_y_continuous(breaks = breaks_pretty(n = 3)) +
         theme_classic()))

(examp.mort <- fullmort %>% 
  filter(thres == 100,
         winter_jday %in% c(168, 217),
         sensor %in% c("ld_caseta_semioberta", "LD_tossal")) %>% 
  mutate(sensor = case_when(sensor == "LD_tossal" ~ "s9 - Ld - O",
                            TRUE ~ "s10 - Ld - O"),
         sensor = factor(sensor,
                         levels = c("s9 - Ld - O", "s10 - Ld - O"))) %>% 
  split(.$winter_jday) %>% 
  map(~ ggplot(data = .x, aes(x = sp_comp, y = mort/100)) +
        geom_point(size = 0.25, alpha = 0.5,
                   position = position_jitterdodge(),
                   aes(color = sp_comp)) +
        geom_boxplot(aes(fill = sp_comp, linetype = sensor),
                     outlier.shape = NA, fatten = 1) +
        scale_color_manual(aesthetics = c("fill", "color"),
                           values = c("P. napi" = "deepskyblue",
                                      "P. rapae" = "goldenrod1"),
                           name = "Species") +
        guides(linetype = "none") +
        facet_wrap(facets = vars(sensor), ncol = 2, scales = "free") +
        labs(x = "Species",
             y = "Mortality")) %>% 
  map2(c(3, 2),
       ~ .x +
         scale_y_continuous(breaks = breaks_pretty(n = .y)) +
         theme_classic() +
         theme(legend.text = element_text(face = "italic"),
               strip.background = element_blank(),
               axis.text.x = element_text(face = "italic"))))

(examp.plot <- examp.temp$`168` + labs(tag = "A") + examp.mort$`168`  +
  examp.temp$`217` + labs(tag = "B") + examp.mort$`217` +
  plot_layout(guides = "collect", ncol = 2))



# Fig. S15: field mortality -----------------------------------------------
sp <- c("P. napi", "P. rapae")
labs <- paste0("<i>", sp, "</i>", " ovipositing sites")
names(labs) <- sp

## A: Daily mortality
(dailymort <- fullmort %>%
  group_by(sensor, winter_jday, sp_comp, thres) %>% 
  summarise(meanmort = mean(mort),
            q975 = quantile(mort, 0.975),
            q025 = quantile(mort, 0.025),
            tmax = mean(tmax),
            tmin = mean(tmin)) %>% 
  left_join(sensor_patch, by = c("sensor" = "Sensor_name")) %>%
  filter(thres == 100) %>% 
  mutate(patch3 = case_when(microhabitat == "O" ~ "O",
                            microhabitat %in% c("SO", "SC") ~ "OC",
                            TRUE ~ NA_character_)) %>% 
  filter((patch3 == "O" & sp_comp == "P. rapae") |
           (patch3 == "OC" & sp_comp == "P. napi")) %>% 
  group_by(patch3, sp_comp, winter_jday, thres) %>% 
  summarise(meanmort = mean(meanmort, na.rm = T)/100,
            tmax = max(tmax, na.rm = T),
            tmin = min(tmin, na.rm = T)) %>% 
  ggplot(aes(x = winter_jday, y = meanmort)) +
  geom_ribbon(aes(ymax = tmax/50, ymin = tmin/50),
              alpha = 0.4) +
  geom_col(aes(color = sp_comp, fill = sp_comp)) +
  geom_hline(aes(yintercept = 0.824), linetype = "dashed") +
  labs(y = "Daily mortality",
       x = "Julian day") +
  facet_wrap(facets = vars(sp_comp), ncol = 2,
             labeller = labeller(sp_comp = labs)) +
  scale_x_continuous(breaks = seq(from = 60, to = 240, by = 60)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.2),
                     sec.axis = sec_axis(trans = ~ .*50,
                                         name = "Temperature (ºC)")) +
  scale_color_manual(aesthetics = c("fill", "color"),
                     values = c("P. napi" = "deepskyblue",
                                "P. rapae" = "goldenrod1"),
                     name = "Species") +
  theme_classic() +
  theme(strip.text = element_markdown(),
        legend.text = element_text(face = "italic")))

## B: mortality during development with TAB
(dvtmort <- fullmort %>% 
  group_by(sensor, sp_comp, thres, simulation) %>% 
  arrange(winter_jday, .by_group = T) %>% 
  mutate(prob_surv = surv/100,
         log_prob = log10(prob_surv),
         log_prob = if_else(log_prob == -Inf, -999, log_prob),
         log_prob_30 = rollsum(log_prob, k = 30, fill = NA, align = "left"),
         prob_30 = 10^log_prob_30) %>% 
  filter(!is.na(prob_30)) %>% 
  group_by(sensor, sp_comp, winter_jday, thres) %>% 
  summarise(mean_prob_30 = mean(prob_30),
            upp_prob_30 = quantile(prob_30, 0.975),
            low_prob_30 = quantile(prob_30, 0.025)) %>% 
  left_join(sensor_patch, by = c("sensor" = "Sensor_name")) %>% 
  mutate(patch3 = case_when(microhabitat == "O" ~ "O",
                            microhabitat %in% c("SO", "SC") ~ "OC",
                            TRUE ~ NA_character_)) %>% 
  filter((patch3 == "O" & sp_comp == "P. rapae") |
           (patch3 == "OC" & sp_comp == "P. napi")) %>% 
  group_by(patch3, sp_comp, winter_jday, thres) %>% 
  summarise(across(.cols = c(contains("prob_30")), mean)) %>% 
  ggplot(aes(x = winter_jday, y = mean_prob_30,
             color = as.factor(thres), fill = as.factor(thres))) +
  geom_line() +
  geom_ribbon(aes(ymax = upp_prob_30, ymin = low_prob_30),
              alpha = 0.25,
              linetype = 0) +
  scale_x_continuous(breaks = seq(from = 60, to = 240, by = 60)) +
  scale_y_continuous(limits = c(0,1)) +
  facet_grid(cols = vars(sp_comp), labeller = labeller(sp_comp = labs)) +
  guides(color = guide_legend(title = "TAB\nthreshold\n(ºC)",
                              title.hjust = 0.5,
                              override.aes =list(fill = NA)),
         fill = "none") +
  labs(x = "Julian day",
       y = "Survival probability\nduring development") +
  theme_classic() +
  theme(strip.text = element_markdown()))

## C: interspecific ration of mean mortality in the field
order_sens <- fullmort %>% 
  group_by(sensor, sp_comp, thres, simulation) %>% 
  arrange(winter_jday, .by_group = T) %>% 
  mutate(prob_surv = surv/100,
         log_prob = log10(prob_surv),
         log_prob = if_else(log_prob == -Inf, -999, log_prob),
         log_prob_30 = rollsum(log_prob, k = 30, fill = NA, align = "left"),
         prob_30 = 10^log_prob_30,
         mort_30 = 1-prob_30,
         high_mort_5 = if_else(mort_30 >= 0.05, 1, 0),
         high_mort_1 = if_else(mort_30 >= 0.01, 1, 0),
         mort_30_filt = if_else(high_mort_1 == 1, mort_30, NA_real_)) %>% 
  filter(!is.na(prob_30),
         sensor != "s8_lepidium") %>% 
  group_by(sensor, sp_comp, thres, simulation) %>% 
  summarise(max_mort_30 = max(mort_30),
            tmax = max(tmax),
            days_highmort_1 = sum(high_mort_1),
            mean_mort_1 = mean(mort_30_filt, na.rm = T),
            mean_mort = mean(mort_30)) %>% 
  filter(thres == 100) %>% 
  group_by(sensor) %>% 
  summarise(mean_mort_1 = mean(mean_mort_1, na.rm = T)) %>% 
  arrange(mean_mort_1) %>% 
  left_join(sensor_patch, by = c("sensor" = "Sensor_name")) %>% 
  filter(!is.na(mean_mort_1), microhabitat != "C", !(microhabitat == "SC" & site == "Me")) %>% 
  rownames_to_column() %>% 
  mutate(site = if_else(site == "AE", "Ld - ", "Me - "),
         site_patch = paste0("s", rowname, " - ",site, microhabitat),
         patch = if_else(microhabitat == "O", "O", "OC"),
         patch = factor(patch, levels = c("O", "OC"))) %>% 
  dplyr::select(-mean_mort_1)

sp <- c("P. napi", "P. rapae")
labs2 <- paste0("<i>", sp, "</i> ovipositing sites")
names(labs2) <- c("OC", "O")

(sensmeanmort.rat <- fullmort %>% 
  group_by(sensor, sp_comp, thres, simulation) %>% 
  arrange(winter_jday, .by_group = T) %>% 
  mutate(prob_surv = surv/100,
         log_prob = log10(prob_surv),
         log_prob = if_else(log_prob == -Inf, -999, log_prob),
         log_prob_30 = rollsum(log_prob, k = 30, fill = NA, align = "left"),
         prob_30 = 10^log_prob_30,
         mort_30 = 1-prob_30,
         high_mort_5 = if_else(mort_30 >= 0.05, 1, 0),
         high_mort_1 = if_else(mort_30 >= 0.01, 1, 0),
         mort_30_filt = if_else(high_mort_1 == 1, mort_30, NA_real_)) %>% 
  filter(!is.na(prob_30),
         sensor != "s8_lepidium") %>% 
  group_by(sensor, sp_comp, thres, simulation) %>% 
  summarise(max_mort_30 = max(mort_30),
            tmax = max(tmax),
            days_highmort_1 = sum(high_mort_1),
            mean_mort_1 = mean(mort_30_filt, na.rm = T),
            mean_mort = mean(mort_30)) %>% 
  filter(thres == 100) %>% 
  pivot_wider(id_cols = c(sensor, thres, simulation),
              names_from = sp_comp,
              values_from = mean_mort_1) %>% 
  rename(PN = "P. napi", PR = "P. rapae") %>% 
  mutate(rat_mort = PN/PR,
         dif_mort = PN - PR) %>% 
  left_join(order_sens) %>% 
  filter(!is.na(site_patch)) %>% 
  mutate(sensor = factor(sensor,
                         levels = order_sens$sensor),
         site_patch = factor(site_patch,
                             levels = order_sens$site_patch),
         patch = factor(patch, levels = c("OC", "O"))) %>% 
  ggplot(aes(x = site_patch, y = rat_mort)) +
  geom_point(size = 0.25, alpha = 0.5,
             position = position_jitter()) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(vars(patch), labeller = labeller(patch = labs2), scales = "free_x") +
  scale_x_discrete(drop = T) +
  labs(y = "Ratio of mean<br>thermal mortality<br>(Pn/Pr)",
       x = "Sensor - site - microhabitat") +
  geom_hline(aes(yintercept = 1)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.title.y = element_markdown(),
        legend.text = element_text(face = "italic"),
        strip.text = element_markdown(),
        strip.background = element_blank()))


(fieldmort.plot <- dailymort + dvtmort + sensmeanmort.rat +
  plot_layout(ncol = 1,
              guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(strip.background = element_blank()))

