############
## Fig. 5 ##
############


# Packages ----------------------------------------------------------------
library(tidyverse)
library(ggtext)
library(patchwork)
library(broom)
library(scales)
library(ggnewscale)
library(zoo)
library(ggdist)


# Data --------------------------------------------------------------------
## Data from TDT experiments
realtdtdata <- read.table("data/TDT_experiment.txt", header = T,
                          sep = ",", dec = ",")

## Data of simulated thermal mortality
fullmort <- read_csv("data/sim_mortalities/mort_min_tdtsim_sensors_full2021-06-25.csv") %>% 
  bind_rows(read_csv("data/sim_mortalities/mort_min_tdtsim_sensors_full2021-06-24.csv")) %>% 
  bind_rows(read_csv("data/sim_mortalities/mort_min_tdtsim_sensors_full2021-06-23.csv"))

## Mean annual series of microclimatic field records at min resolution
sensors_min_spline <- readRDS("data/sim_mortalities/sensors_minute_spline_2021-06-25.RDS")

dayjday <- fullmort %>% 
  distinct(sensor, winter_jday) %>% 
  group_by(sensor) %>% 
  mutate(day = min_rank(winter_jday))

sensors_spline <- sensors_min_spline %>% 
  mutate(day = floor((hour-1)/24) + 1) %>% 
  left_join(dayjday)
rm(sensors_min_spline)

# correspondence of each sensor to a microhabita type
sensor_patch <- read_delim("data/table_yearsens.txt",
                                           delim = " ",
                                           col_names = T)



# A: tdt curve ------------------------------------------------------------
tdtinters <- realtdtdata %>% 
  nest(data = -sp) %>% 
  mutate(lm = map(data,
                  ~lm(log10(aprox_minute_dead) ~ SENSOR_mean_temp,
                      data = .)),
         est = map(lm, tidy),
         est = map(est, dplyr::select, term, estimate)) %>% 
  unnest(est) %>% 
  dplyr::select(sp, term, estimate) %>% 
  mutate(term = c("int", "slope", "int", "slope")) %>% 
  pivot_wider(names_from = c(sp, term),
              names_sep = "_",
              values_from = estimate) %>% 
  mutate(temp = (PN_int - PR_int)/(PR_slope - PN_slope),
         time = 10^(PN_int + PN_slope*temp)) %>% 
  add_row(temp = rep(.$temp, 249),
          time = rep(.$time, 249)) %>% 
  mutate(temp_y = 7:256,
         time_x = seq(39, 41.2, length.out = 250))

(tdtplot <- realtdtdata %>% 
  ggplot(aes(x = SENSOR_mean_temp, y = aprox_minute_dead,
             color = sp, fill = sp)) +
  scale_y_log10(breaks = c(10, 60, 180, 420, 900),
                expand = expansion(mult = c(0, 0))) +
  scale_color_manual(aesthetics = c("fill", "color"),
                     values = c("PN" = "deepskyblue",
                                "PR" = "goldenrod1"),
                     name = "Species",
                     labels = c("P. napi", "P. rapae")) +
  scale_x_continuous(breaks = c(40, 42, 44),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(y = "Survival time (min, log scale)",
       x = "Temperature (ºC)") +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  geom_line(data = tdtinters, aes(x = temp, y = temp_y,
                                  color = NULL, fill = NULL),
            linetype = "dashed") +
  geom_line(data = tdtinters, aes(x = time_x, y = time,
                                  color = NULL, fill = NULL),
            linetype = "dashed") +
  geom_text(aes(x = 41, y = 10, label = "T*", xjust = 1, yjust = 1, size = 8,
                fontface = "bold.italic"), color = "black") +
  theme_classic() +
  theme(legend.text = element_text(face = "italic"),
        axis.title.x = element_text(size = 10)))


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
           mort_30_filt = if_else(high_mort_1 == 1, mort_30, NA_real_)) %>% 
    filter(!is.na(prob_30),
           sensor != "s8_lepidium") %>% 
    group_by(sensor, sp_comp, thres, simulation) %>% 
    summarise(max_mort_30 = max(mort_30),
              tmax = max(tmax),
              days_highmort_1 = sum(high_mort_1),
              mean_mort_1 = mean(mort_30_filt, na.rm = T),
              mean_mort = mean(mort_30)) %>% 
    left_join(order_sens) %>% 
    filter(!is.na(site_patch)) %>% 
    mutate(sensor = factor(sensor,
                           levels = order_sens$sensor),
           site_patch = factor(site_patch,
                               levels = order_sens$site_patch),
           patch = factor(patch, levels = c("OC", "O")),
           side_sp = if_else(sp_comp == "P. napi", "left", "right")) %>% 
    ggplot(aes(x = site_patch, y = mean_mort_1)) +
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



# C: interspecific ratio of mean mortality in the field -------------------
(sensmeanmort.rat <- fullmort %>% 
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
    # geom_point(size = 0.25, alpha = 0.5,
    #            position = position_jitter()) +
   # geom_violin(scale = "width") +
   stat_eye(
     #aes(fill = stat(y > 1)),
     fill = "grey",
     normalize = "groups",
     adjust = .5,
     width = .75,
     alpha = .75,
     .width = 0,
     point_colour = NA) +
   geom_boxplot(outlier.shape = NA,
                width = 0.25,
                 fill = NA,
                coef = 0) +
    # facet_wrap(vars(patch), labeller = labeller(patch = labs2), scales = "free_x") +
   facet_wrap(vars(patch),
              labeller = labeller(patch = labs2),
              scales = "free") + 
   scale_x_discrete(drop = T) +
   scale_y_continuous(breaks = seq(0.5, 2, by = 0.5),
                      limits = c(0.15, 2.2)) +
    labs(y = "Ratio of mean<br>thermal mortality<br>(*P. napi*/*P. rapae*)",
         x = "Sensor - site - microhabitat") +
    geom_hline(aes(yintercept = 1)) +
   scale_fill_manual(values = c("goldenrod1", "deepskyblue")) +
   guides(fill = "none") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          axis.title.y = element_markdown(),
          axis.title.x = element_text(size = 10),
          legend.text = element_text(face = "italic"),
          strip.text = element_markdown(),
          strip.background = element_blank()))



# assemblage --------------------------------------------------------------
(tdt.mort.plot <- tdtplot +
  labs(y = "Survival time<br>(min, log scale)",
       tag = "A") +
  guides(color = "none", fill = "none", size = "none") +
  sensmeanmort +
  labs(y = "Mean thermal<br>mortality",
       tag = "B") +
   theme(axis.title.x = element_blank(),
         axis.text.x = element_blank()) +
  guides(color = "none",
         fill = guide_legend(override.aes = (list(color = NA)))) +
   sensmeanmort.rat +
   labs(tag = "C") +
   theme(strip.text = element_blank()) +
  plot_layout(ncol = 1, guides = "collect") &
  theme(text = element_text(size = 10),
        plot.background = element_rect(fill = NA),
        plot.tag = element_text(face = "bold", size = 12),
        plot.title = element_text(size = 10),
        axis.title.y = element_markdown(size = 10),
        # axis.title.x = element_text(size = 10),
        # strip.text = element_markdown(size = 10),
        legend.title = element_markdown(size = 10),
        legend.position = "top",
        axis.text.y = element_text(size = 8)))



# ggsave(filename = "figures/fig5.png",
#        plot = tdt.mort.plot,
#        height = 21,
#        width = 15,
#        units = "cm",
#        dpi = 600)
