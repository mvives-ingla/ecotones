######################
# oviposition models #
######################


# packages ----------------------------------------------------------------
library (tidyverse)
library (lme4)
library(performance)


# data --------------------------------------------------------------------
## number of ovipositions per female
ovi.byfem <- read.csv ("data/ovi_by_female.csv")
## number of ovipositing females per survey
ovi.bysurv <- read.csv ("data/ovi_by_survey.csv")


# model 1: number of oviposition per female -------------------------------
# ovi.byfem <- ovi.byfem %>% 
#   rownames_to_column() %>% 
#   mutate(fem_code = if_else(is.na(fem_code), as.integer(rowname), fem_code))

mod1 <- glmer (n_ovi ~
                 microhabitat*species +
                 (1 | day_period) +
                 (1 | jday_cat) +
                 (1 | site) +
                 # (1 | fem_code) +
                 offset (log (duration)), # same as offset = log(duration)
               # na.action = na.pass,
               data = ovi.byfem,
               family = poisson)

summary (mod1)

## extraction of fixed effects coefficients and their statistical significance as a dataframe
fixef1 <- round (coef (summary (mod1), complete = T), digits = 4)

int <- fixef1[1,1] # coefficient for P. napi ovipositions in C microhabitats
efsize <- exp(c(int, fixef1[2:6,1] + int)) # all coefficients at the response scale
percent.pn <- efsize[1:3]/sum(efsize[1:3]) # distribution of ovipositions of P. napi
ratio.pn <- efsize[3]/efsize[1:2] # ratio of ovipositions in OC microhabitats vs C and O microhabitats
efsize.pr <- efsize[4:6]/sum(efsize[4:6]) # distribution of ovipositions of P. napi



## extraction of random effects variance estimate as a dataframe
ranef1 <- as.data.frame (VarCorr (mod1, comp = c ("Variance", "Std.Dev.")))

## test of effects
### Analysis of deviance, Wald chi-sq test
eftest1 <- car::Anova (mod1, type = "III")


r2(mod1)

# Fig. S4 -----------------------------------------------------------------
(ran.plot1 <- ranef (mod1, condVar = T, whichel = c("day_period")) %>% 
  as.data.frame() %>% 
  mutate (grp = factor (grp, levels = c("1", "2", "3", "4", "5"))) %>% 
  ggplot (aes (y = grp, x = exp (condval))) +
  geom_point() +
  geom_errorbarh (aes (xmin = exp (condval - 1.96*condsd),
                       xmax = exp (condval + 1.96*condsd),
                       height = 0)) +
  labs (y = "Period of the day",
        x = "Predicted num. of ovipostions",
        title = "Both sites") +
  scale_y_discrete (labels = c ("<11", "11–13", "13–15", "15–17", ">17")) +
  theme_classic () +
  theme (text = element_text (size = 10),
         axis.title = element_text (size = 10),
         axis.text = element_text (size = 10),
         plot.title = element_text (size = 10,
                                    face = "plain"),
         panel.border = element_rect (fill = NA),
         plot.tag = element_text (size = 10,
                                  face = "bold")) +
  guides (size = "none"))



# model 2: number of ovipositing females by survey ------------------------
mod2 <- glmer (n_fem ~
                 microhabitat*species +
                 (1 | day_period) +
                 (1 | jday_cat) +
                 (1 | site) +
                 offset (log(duration)),
               data = ovi.bysurv,
               family = poisson)

summary (mod2)

## extraction of fixed effects coefficients and their statistical significance as a dataframe
fixef2 <- round (coef (summary (mod2)), digits = 4)

## extraction of random effects variance estimate as a dataframe
ranef2 <- as.data.frame (VarCorr (mod2, comp = c ("Variance", "Std.Dev.")))

## extraction of fixef effects significance test as a dataframe
eftest2 <- car::Anova (mod2, type = "III")

r2(mod2)



