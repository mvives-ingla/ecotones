##############
# TDT models #
##############

# Packages ----------------------------------------------------------------
library (tidyverse)
library(broom)
library(performance)
library(lme4)
library(ggtext)

# Functions ---------------------------------------------------------------
source("code/thermal_mortality/Thermal_landscape_functions_mod.R")


# Data --------------------------------------------------------------------
## Data from TDT experiments
realtdtdata <- read.table("data/TDT_experiment.txt", header = T,
                          sep = ",", dec = ",")


# Models ------------------------------------------------------------------

## ANCOVA: species x treatment to test whether the two species present signficantly different TDT curves
ancova <- lm(log10(aprox_minute_dead) ~ sp*SENSOR_mean_temp, data = realtdtdata) %>% 
  tidy()

ancova$estimate[3]/(ancova$estimate[3] + ancova$estimate[4]) ## slope of the TDT curves of P. napi vs P. rapae

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
         coefs = map(summ, coef),
         confint = map(model, confint),
         r2 = map(model, r2),
         eff = map(model, car::Anova,type = "III"))


mods$summ %>% 
  set_names(nm = mods$sp)

mods$eff %>% 
  set_names(nm = mods$sp)


mods$r2 %>% 
  set_names(nm = mods$sp)

realtdtdata %>% 
  split(.$sp) %>% 
  map(~ ggplot(data = .,
               aes(x = initial.weight, y = log10(aprox_minute_dead),
             color = as.factor(treatment))) +
        geom_point() +
        geom_smooth(method = "lm") +
        geom_smooth(aes(color = NULL), method = "lm") +
    labs(title = .$sp) +
  theme_bw())


realtdtdata %>% 
  split(.$sp) %>% 
  map(~ ggplot(data = .,
               aes(x = SENSOR_mean_temp, y = log10(aprox_minute_dead),
                   color = as.factor(Site))) +
        geom_point() +
        geom_smooth(method = "lm") +
        # geom_smooth(aes(color = NULL), method = "lm") +
        labs(title = .$sp) +
        theme_bw())


### Quantifying the effect of weight

#### The increase in weight needed to increase survival time an order of magnitude
order_mag <- mods %>% 
  mutate(confint = map(confint,
                       as.data.frame),
         coefs = map(coefs,
                     as.data.frame),
         confint = map(confint,
                       ~ .[3:6,]),
         coefs = map(coefs,
                     ~ .[,1])) %>% 
  select(sp, coefs, confint) %>% 
  unnest(c(coefs, confint)) %>% 
  rownames_to_column(var = "term") %>% 
  mutate(term = rep(c("int", "temperature", "weight", "site"), times = 2)) %>% 
  rename(conflow = "2.5 %",
         confhigh = "97.5 %") %>% 
  pivot_longer(-c(term, sp),
               names_to = "estimate",
               values_to = "value") %>% 
  mutate(inverse = 1/value) %>% 
  pivot_wider(id_cols = c(term, sp),
              values_from = inverse,
              names_from = estimate) %>% 
  select(sp, term, conflow, coefs, confhigh)

#### The increase in time of survival when weight increases 0,1g
incr_surv <- mods %>% 
  mutate(confint = map(confint,
                       as.data.frame),
         coefs = map(coefs,
                     as.data.frame),
         confint = map(confint,
                       ~ .[3:6,]),
         coefs = map(coefs,
                     ~ .[,1])) %>% 
  select(sp, coefs, confint) %>% 
  unnest(c(coefs, confint)) %>% 
  rownames_to_column(var = "term") %>% 
  mutate(term = rep(c("int", "temperature", "weight", "site"), times = 2),
         incr_pred = rep(c(0, 1, 0.1, 1), times = 2)) %>% 
  rename(conflow = "2.5 %",
         confhigh = "97.5 %") %>% 
  pivot_longer(-c(term, sp, incr_pred),
               names_to = "estimate",
               values_to = "value") %>% 
  mutate(rat_t2t1 = (10^(value*incr_pred) - 1)*100) %>% 
  pivot_wider(id_cols = c(term, sp),
              values_from = rat_t2t1,
              names_from = estimate) %>% 
  select(sp, term, conflow, coefs, confhigh) %>% 
  mutate(incr = case_when(term == "temperature" ~ "+1ºC<br>in temperature",
                          term == "weight" ~ "+0.1g<br>in larval weight",
                          term == "site" ~ "from lowland<br>to mid-elevation",
                          T ~ NA_character_))



(plot.effects <- incr_surv %>% 
  filter(!is.na(incr)) %>% 
  ggplot(aes(x = incr, y = coefs)) +
  geom_col(aes(fill = sp), alpha = .5, position = position_dodge(width = 1)) +
  geom_linerange(aes(ymin = conflow, ymax = confhigh, color = sp),
                 position = position_dodge(width = 1)) +
  labs(x = "Increase in predictor",
       y = "Relative change in survival time (%)") +
  scale_fill_manual(values = c("deepskyblue", "goldenrod1"),
                    aesthetics = c("fill", "color"),
                    name = "Species",
                    labels = c("<i>Pieris napi</i>", "<i>Pieris rapae</i>")) +
  # scale_y_continuous(labels = c("-50", "+0", "+50", "+100")) +
  theme_classic() +
  theme(axis.text.x = element_markdown(),
        legend.text = element_markdown()))


