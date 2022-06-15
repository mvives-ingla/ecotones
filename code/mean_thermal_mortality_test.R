library(tidyverse)
library(zoo)


## Data of simulated thermal mortality
fullmort <- read_csv("data/sim_mortalities/mort_min_tdtsim_sensors_full2021-06-25.csv") %>% 
  bind_rows(read_csv("data/sim_mortalities/mort_min_tdtsim_sensors_full2021-06-24.csv")) %>% 
  bind_rows(read_csv("data/sim_mortalities/mort_min_tdtsim_sensors_full2021-06-23.csv"))

# correspondence of each sensor to a microhabita type
sensor_patch <- read_delim("data/table_yearsens.txt",
                           delim = " ",
                           col_names = T)


fullmort %>% 
  filter(thres == 100) %>% 
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
         mort_30_filt = if_else(high_mort_1 == 1, mort_30, NA_real_),
         mort_30_summ = if_else(winter_jday %in% 150:245, mort_30, NA_real_)) %>% 
  filter(!is.na(prob_30),
         sensor != "s8_lepidium") %>% 
  group_by(sensor, sp_comp, thres, simulation) %>% 
  summarise(max_mort_30 = max(mort_30),
            tmax = max(tmax),
            days_highmort_1 = sum(high_mort_1),
            mean_mort_1 = mean(mort_30_filt, na.rm = T),
            mean_mort = mean(mort_30),
            mean_mort_summ = mean(mort_30_summ, na.rm = T)) %>% 
  pivot_longer(cols = mean_mort_summ:mean_mort_1,
               names_to = "type",
               values_to = "mort") %>% 
  split(.$sensor) %>% 
  map(~ ggplot(data = .,
               aes(x = sp_comp, y = mort, fill = type)) +
  geom_violin() +
  labs(title = .$sensor) +
  theme_classic())
## La mortalitat mitjana agafant tots els valors de l'estiu i l'anual només agafant les acumulades >1%
# son generalment similars, tot i que no sempre.
## A vegades mortalitats d'estiu més elevades que les mortalitats > 1%, especialment en els llocs oberts. Probablement pq els períodes amb mortalitat > 1% van més enllà de l'estiu
## En llocs semi-oberts, com que els períodes amb mortalitat es concentren a l'estiu, les diferències son menors
## Altres variacions poden ser degudes a la diferent inclusió de zeros en els càlculs



fullmort %>% 
  filter(thres == 100) %>% 
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
         mort_30_filt = if_else(high_mort_1 == 1, mort_30, NA_real_),
         mort_30_summ = if_else(winter_jday %in% 150:245, mort_30, NA_real_)) %>% 
  filter(!is.na(prob_30),
         sensor != "s8_lepidium") %>% 
  # group_by(sensor, sp_comp, thres, simulation) %>% 
  # summarise(max_mort_30 = max(mort_30),
  #           tmax = max(tmax),
  #           days_highmort_1 = sum(high_mort_1),
  #           mean_mort_1 = mean(mort_30_filt, na.rm = T),
  #           mean_mort = mean(mort_30),
  #           mean_mort_summ = mean(mort_30_summ, na.rm = T)) %>% 
  select(c(sensor:sp_comp, thres, simulation, mort_30, mort_30_filt, mort_30_summ)) %>% 
  pivot_longer(cols = mort_30:mort_30_summ,
               names_to = "type",
               values_to = "mort") %>% 
  split(.$sensor) %>% 
  map(~ ggplot(data = .,
               aes(x = mort, color = type)) +
        geom_density(aes(x = mort, after_stat(scaled))) +
        labs(title = .$sensor) +
        facet_wrap(vars(sp_comp), nrow = 2, scales = "free_y") +
        theme_classic())



# B: mean mortality in the field ------------------------------------------
order_sens <- fullmort %>% 
  filter(thres == 100) %>% 
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
  group_by(sensor) %>% 
  summarise(mean_mort_1 = mean(mean_mort_1, na.rm = T)) %>% 
  #mutate(mean_mort_5 = if_else(is.na(mean_mort_5), 0, mean_mort_5)) %>% 
  arrange(mean_mort_1) %>% 
  left_join(sensor_patch, by = c("sensor" = "Sensor_name")) %>% 
  filter(!is.na(mean_mort_1), patch2 != "C", !(patch2 == "SC" & site == "Me")) %>% 
  rownames_to_column() %>% 
  mutate(site = if_else(site == "Ld", "Ld - ", "Me - "),
         site_patch = paste0("s", rowname, " - ",site, patch2),
         patch = if_else(patch2 == "O", "O", "OC"),
         patch = factor(patch, levels = c("O", "OC"))) %>% 
  dplyr::select(-mean_mort_1)

sp <- c("P. napi", "P. rapae")
labs2 <- paste0("Microhabitat selected by <i>", sp, "</i>")
names(labs2) <- c("OC", "O")


(sensmeanmort <- fullmort %>% 
    filter(thres == 100) %>% 
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
           mort_30_filt = if_else(high_mort_1 == 1, mort_30, NA_real_),
           mort_30_summ = if_else(winter_jday %in% 150:245, mort_30, NA_real_)) %>% 
    filter(!is.na(prob_30),
           sensor != "s8_lepidium") %>% 
    group_by(sensor, sp_comp, thres, simulation) %>% 
    summarise(max_mort_30 = max(mort_30),
              tmax = max(tmax),
              days_highmort_1 = sum(high_mort_1),
              mean_mort_1 = mean(mort_30_filt, na.rm = T),
              mean_mort = mean(mort_30),
              mean_mort_summ = mean(mort_30_summ, na.rm = T)) %>% 
    left_join(order_sens) %>% 
    filter(!is.na(site_patch)) %>% 
    mutate(sensor = factor(sensor,
                           levels = order_sens$sensor),
           site_patch = factor(site_patch,
                               levels = order_sens$site_patch),
           patch = factor(patch, levels = c("OC", "O")),
           side_sp = if_else(sp_comp == "P. napi", "left", "right")) %>% 
    ggplot(aes(x = site_patch, y = mean_mort_summ)) +
    ggdist::stat_halfeye(
      aes(side = side_sp, fill = sp_comp, color = sp_comp),
      normalize = "groups",
      alpha = .75,
      adjust = .5,
      width = .5,
      # treure l'interval
      .width = 0,
      point_colour = NA
    ) +
    geom_boxplot(aes(group = paste(site_patch, sp_comp)),
                 size = 0.4,
                 width = .25,
                 fill = NA,
                 outlier.shape = NA,
                 coef = 0) +
    facet_wrap(vars(patch), labeller = labeller(patch = labs2), scales = "free") +
    scale_x_discrete(drop = T) +
    scale_color_manual(aesthetics = c("fill", "color"),
                       values = c("P. napi" = "deepskyblue",
                                  "P. rapae" = "goldenrod1"),
                       name = "Species") +
    labs(y = "Mean thermal mortality",
         x = "Sensor - site - microhabitat") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.text = element_text(face = "italic"),
          strip.text = element_markdown(size = 10),
          strip.background = element_blank()))
