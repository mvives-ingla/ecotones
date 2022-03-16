#Data preparation for Proceedings B

library (tidyverse)
library (readxl)
library (lubridate)
library (hms)

#############################################
# Microenvironment and host plant variables #
#############################################

## Basic dataset

data.leaf <- read_excel ("C:/Users/MariaCyborg/OneDrive - CREAF/Dades i analisis/Bases definitives/Host_plant_2017/FIELD_DATASET_2017_v8.xlsx", sheet = 2)
data.ind <- read_excel ("C:/Users/MariaCyborg/OneDrive - CREAF/Dades i analisis/Bases definitives/Host_plant_2017/FIELD_DATASET_2017_v8.xlsx", sheet = 3)
data.patch <- read_excel ("C:/Users/MariaCyborg/OneDrive - CREAF/Dades i analisis/Bases definitives/Host_plant_2017/FIELD_DATASET_2017_v8.xlsx", sheet = 4)


data <- list (
  leaf = data.leaf,
  ind = data.ind, 
  patch = data.patch
)


data <- data %>% 
  map (~ mutate_if (., is.character, as.factor)) %>% 
  map (~ unite (., "cond_obv_temp",
                c("rad_obv_temp", "wind_obv_temp"), remove = F)) %>% 
  map (~ mutate (.,
                 patch_plot = recode (patch_plot,
                                      "O - LD41" = "O - LD4",
                                      "OO - LD42" = "O - LD4"),
                 leaf_dens = .$dry_weight*1000/.$length))

rm (data.patch, data.ind, data.leaf)

# Random subset

set.seed (7)

data.sample <- data [["ind"]] %>%
  filter (seguiment == "fixed") %>%
  group_by (individual) %>% 
  sample_n (1)

data.ran <- data [["ind"]] %>% 
  filter (seguiment == "once") %>%
  bind_rows (data.sample)

data.ran <- list (ind = data.ran)

#data.sample.leaf <- data [["leaf"]] %>%
 # filter (seguiment == "fixed") %>%
 # select (individual, jday) %>%
 # group_by (individual) %>%
 # sample_n (1) %>% #hauria d'haver seleccionat els mateixos individus que pel nivell d'individus, no una altra mostra a l'atzxar
 # left_join (data [["leaf"]], by = c("jday", "individual"))

## si volgués que les bases de dades quadressin entre nivells hauria de ser així
data.sample.leaf <- data.sample %>% 
  dplyr::select(jday, individual) %>% 
  left_join (data [["leaf"]], by = c("jday", "individual"))

identical(distinct(data.sample, individual, jday), distinct(data.sample.leaf, individual, jday))



data.ran.leaf <- data [["leaf"]] %>% 
  filter (seguiment == "once") %>%
  bind_rows (data.sample.leaf)


data.ran$leaf <- data.ran.leaf
data.ran$patch <- data [["patch"]]

rm(list = c("data.sample", "data.sample.leaf", "data.ran.leaf", "data.sample.leaf.2"))



## Data.patch preparation

CJ2.113 <- data.ran[["patch"]] %>% 
  filter (jday == 113,
          patch == "CJ2") %>% 
  dplyr::select (patch, starts_with ("ori"))

CJ2newrows <- data.ran[["patch"]] %>%
  rownames_to_column () %>% 
  mutate (rowname = as.numeric (rowname)) %>% 
  filter (jday > 113,
          jday <= 176,
          patch == "CJ2") %>% 
  dplyr::select (jday, patch) %>% 
  left_join (CJ2.113, by = c("patch" = "patch")) # change needed values in cover (forest CJ)


data.microh <- data[["patch"]] %>% 
  mutate (soil_surface_R = case_when ((jday == 206 | ### removing outliers of surface temperatures (don't remember why are removed: 7_thermal conditions and my notebooks may have the answer)
                                         jday == 167 |
                                         jday == 153 |
                                         jday == 159 |
                                         jday == 164) ~ NA_real_,
                                      TRUE ~ soil_temp_R_fback),
          soil_surface_NR = case_when ((jday == 206 |
                                          jday == 167 |
                                          jday == 153 |
                                          jday == 159 |
                                          jday == 164) ~ NA_real_,
                                       TRUE ~ soil_temp_NR_fback),
          soil_surface = case_when ((!is.na(soil_surface_R) & !is.na(soil_surface_NR)) ~
                                         (soil_surface_R + soil_surface_NR)/2,
                                    (!is.na(soil_surface_R) & is.na(soil_surface_NR)) ~
                                         soil_surface_R,
                                    (is.na(soil_surface_R) & !is.na(soil_surface_NR)) ~ 
                                         soil_surface_NR,
                                    TRUE ~ NA_real_),
          soil_temp_hour = case_when (is.na (soil_surface) ~ NA_real_,
                                      TRUE ~ soil_temp_hour),
          outlier_soil_fb = if_else (is.na (soil_surface), NA_real_, outlier_soil_fb)) %>% # generate soil surface temperatures
  mutate_at (vars(starts_with ("ori"), contains ("cover")),
             ~ case_when ((patch == "LD5" | jday == 70 | jday == 206 | jday == 207) ~ NA_real_, 
                          TRUE ~ .)) %>% # delete LD5 values of cover and outliers
  rows_update (CJ2newrows, by = c("jday", "patch")) %>% # change needed values in cover (forest CJ)
  mutate (microhabitat = recode (patch, "LD1" = "SC", "LD2" = "C", "LD3" = "SO", "LD41" = "O", "LD42" = "O",
                                 "LD5" = "SC", "CJ1" = "SC", "CJ2" = "C", "CJ3" = "SO", "CJ4" = "O",
                                 .default = NA_character_)) %>% # modificar patch i LD5 --> SC
  filter (species != "Brassica nigra",
          !is.na (microhabitat)) %>% 
  dplyr::select (site,jday, month, month_name, season_name, species, microhabitat, #identifiers
          soil_surface_R, soil_surface_NR, soil_surface, #soil temperature
          soil_temp_hour, METCAT_TM, METCAT_TX, SENS_TM, SENS_TX, SENS_TM2H, outlier_soil_fb, #accompanying soil temp
          ground_cover, hei_ground_cover, starts_with ("ori"), #canopy and ground cover
          n_resprouts, #resprout density
          n_veg, n_rep, n_sen, n_sum_ros) %>% #phenology
  rename (resprout_dens = n_resprouts) %>% 
  filter_at (vars(soil_surface_R:soil_temp_hour, outlier_soil_fb:n_sum_ros), any_vars (!is.na (.))) 

write_csv (data.microh, "data/data_microhabitat.csv")


rm (CJ2.113, CJ2newrows)


## Data.ind preparation

ind.codes <- data [["ind"]] %>% 
  mutate (microhabitat = recode (patch, "LD1" = "SC", "LD2" = "C", "LD3" = "SO", "LD41" = "O", "LD42" = "O",
                                 "LD5" = "SC", "CJ1" = "SC", "CJ2" = "C", "CJ3" = "SO", "CJ4" = "O",
                                 .default = NA_character_)) %>% 
  dplyr::select (site, jday, patch, microhabitat, individual) %>% 
  filter (!is.na (microhabitat)) %>% 
  group_by (site, microhabitat) %>% 
  mutate (row = row_number()) %>% 
  unite ("ind", site, microhabitat, row, remove = F) %>% 
  dplyr::select (-row)

data.plants.ran <- data.ran [["ind"]] %>% 
  dplyr::select (jday, patch, individual, phenology, cycle_phase, height, hei.with.leaves, n_leaves)

data.plants <- data [["ind"]] %>%
  mutate (soil_temp_couple = case_when ((jday == 206 | ### removing outliers of soil temperatures (don't remember why are removed: 7_thermal conditions and my notebooks may have the answer)
                                         jday == 167 |
                                         jday == 153 |
                                         jday == 159 |
                                         jday == 164) ~ NA_real_,
                                      TRUE ~ soil_temp_couple),
          air_temp = case_when ((jday == 206 |
                                   jday == 167 |
                                   jday == 153 |
                                   jday == 159 |
                                   jday == 164) ~ NA_real_,
                                TRUE ~ air_mean),
          soil_temp_hour = if_else (is.na (soil_temp_couple), NA_real_, soil_temp_hour),
          air_hour = if_else (is.na (air_temp), NA_real_, air_hour),
          rad_air = if_else (is.na (air_temp), NA_character_, as.character(rad_air)),
          wind_air = if_else (is.na (air_temp), NA_character_, as.character(wind_air)),
          microhabitat = recode (patch, "LD1" = "SC", "LD2" = "C", "LD3" = "SO", "LD41" = "O", "LD42" = "O",
                                 "LD5" = "SC", "CJ1" = "SC", "CJ2" = "C", "CJ3" = "SO", "CJ4" = "O",
                                 .default = NA_character_)) %>% # modificar patch i LD5 --> SC
  dplyr::select (site, jday, month, month_name, season_name, species, microhabitat,
          patch, individual, #identifiers
          soil_temp_couple, soil_temp_hour, air_temp, air_hour, rad_air, wind_air, METCAT_TM, METCAT_TX, SENS_TM, SENS_TX, SENS_TM2H,
          org, min, volts) %>% 
  left_join (data.plants.ran) %>% 
  left_join (ind.codes) %>% 
  dplyr::select (site:microhabitat, ind, soil_temp_couple:n_leaves) %>% 
  rename (stem_length = height,
          soil_hum_org = org,
          soil_hum_min = min,
          soil_hum_volts = volts,
          resp_height = hei.with.leaves,
          resp_leaves = n_leaves) %>% 
  filter (!is.na (microhabitat), species != "Brassica nigra") %>% 
  filter_at (vars(soil_temp_couple:wind_air,
                  soil_hum_org:soil_hum_volts,
                  stem_length:resp_leaves),
             any_vars (!is.na (.))) 
          
write_csv (data.plants, "data/data_plants.csv")

rm (data.plants.ran)

## data.leaves preparation

data.leaves.ran <- data.ran [["leaf"]] %>% 
  mutate (fresh_weight = case_when (is.na(filt_wcnorm) ~ NA_real_,
                                    individual == "CJ2_24" & dry_weight < 0.001 ~ NA_real_,
                                    TRUE ~ fresh_weight),
          dry_weight = case_when (is.na(filt_wcnorm) ~ NA_real_,
                                    individual == "CJ2_24" & dry_weight < 0.001 ~ NA_real_,
                                    TRUE ~ dry_weight),
          water_content = (fresh_weight - dry_weight)/dry_weight,
          leaf_dens = if_else(is.na(water_content), NA_real_, leaf_dens)) %>%
  dplyr::select (jday, patch, individual, length, width, chl, fresh_weight, dry_weight, water_content,
          leaf_dens, leaf_code)


data.leaves <- data [["leaf"]] %>% 
  mutate (obv_temp = case_when ((jday == 206 | ### removing outliers of soil temperatures (don't remember why are removed: 7_thermal conditions and my notebooks may have the answer)
                                           jday == 167 |
                                           jday == 153 |
                                           jday == 159 |
                                           jday == 164) ~ NA_real_,
                                        TRUE ~ obv_temp),
          rev_temp = case_when ((jday == 206 |
                                   jday == 167 |
                                   jday == 153 |
                                   jday == 159 |
                                   jday == 164) ~ NA_real_,
                                TRUE ~ rev_temp),
          rad_obv_temp = if_else (is.na (obv_temp), NA_character_, as.character(rad_obv_temp)),
          wind_obv_temp = if_else (is.na (obv_temp), NA_character_, as.character(wind_obv_temp)),
          cond_obv_temp = if_else (is.na (obv_temp), NA_character_, as.character (cond_obv_temp)),
          microhabitat = recode (patch, "LD1" = "SC", "LD2" = "C", "LD3" = "SO", "LD41" = "O", "LD42" = "O",
                                 "LD5" = "SC", "CJ1" = "SC", "CJ2" = "C", "CJ3" = "SO", "CJ4" = "O",
                                 .default = NA_character_)) %>% 
  dplyr::select (site, jday, month, month_name, season_name, species, microhabitat,
          patch, individual, phenology, cycle_phase, leaf_type, leaf_code, #identifiers)
          leaf_state, obv_temp:wind_obv_temp, METCAT_TM, METCAT_TX, SENS_TM, SENS_TX, SENS_TM2H) %>%
  left_join (data.leaves.ran) %>% 
  left_join (ind.codes) %>% 
  dplyr::select (site:microhabitat, ind, phenology:leaf_dens) %>% 
  filter (!is.na (microhabitat), species != "Brassica nigra") %>% 
  filter_at (vars (obv_temp:rev_temp, length:leaf_dens), any_vars (!is.na (.)))


write_csv (data.leaves, "data/data_leaves.csv")

rm (data.leaves.ran)


## butterfly phenology

but.pheno <- read_excel ("C:/Users/MariaCyborg/OneDrive - CREAF/Dades i analisis/cBMS/arxius_originals/napi_rapae_2017-2018_zeros.xls", sheet = 1) %>% 
  mutate (jday = yday (Datam)) %>% 
  filter (Any == "2017") %>% 
  mutate (Ab_index = case_when (IDitin == "09" ~ SumaDeNindiv/1672*1000,
                                IDitin == "01" ~ SumaDeNindiv/4296*1000)) %>% 
  rename (site = IDitin) %>% 
  select (site, IDesp, jday, Ab_index) %>% 
  mutate (site = recode (site, `01` = "AE", `09` = "CJ"))


write_csv (but.pheno, "final_versions/ecotones_procB/but_pheno.csv")


## Fruits dataset


ld <- read_excel ("C:/Users/MariaCyborg/OneDrive - CREAF/Dades i analisis/Bases definitives/Host_plant_2017/Fruits_DEF_v3_2019.02.19.xlsx", sheet = 1) %>% 
  mutate (microhabitat = recode (patch, "LD1" = "SC", "LD2" = "C", "LD3" = "SO", "LD41" = "O", "LD42" = "O",
                                 .default = NA_character_),
          species = "Lepidium draba") %>% 
  rename (n_fruits = total_fruits) %>% 
  select (species, microhabitat, n_fruits)


ap <- read_excel ("C:/Users/MariaCyborg/OneDrive - CREAF/Dades i analisis/Bases definitives/Host_plant_2017/Fruits_DEF_v3_2019.02.19.xlsx", sheet = 3) %>% 
  mutate (microhabitat = recode (patch, "CJ1" = "SC", "CJ2" = "C", "CJ3" = "SO", "CJ4" = "O",
                                 .default = NA_character_),
          species = "Alliaria petiolata") %>% 
  select (species, microhabitat, n_fruits)

fruits.data <- bind_rows (ap, ld) %>% 
  filter (!is.na (microhabitat))

write_csv (fruits.data, "final_versions/ecotones_procB/data_fruits.csv")


###############
# Oviposition #
###############

## oviposition/female: model 6.4 (main text)
############################################

#per a separar els ous ovipositats per cada femella (efecte femella) i, a la vegada, tenir zeros.
fems <- read_excel ("C:/Users/MariaCyborg/OneDrive - CREAF/Dades i analisis/Bases definitives/Behaviour_oviposition/Behaviour_20172018_20190723.xlsx",
                    sheet = 2) %>% 
  mutate (Dia = as.numeric (Dia),
          Dia = as_date (Dia, origin = "1899-12-30")) %>% 
  rename (h_ini = "H inicial",
          h_fin = "H final",
          species = "Espècie papallona",
          individual = "Número individu",
          plot = "Parcel·la") %>% 
  mutate (h_ini = as_hms (h_ini),
          h_fin = as.numeric (h_fin),
          h_fin = h_fin*60*60*24,
          h_fin = as_hms (h_fin)) %>% 
  filter (Month >= 6,
          Month <= 9) %>% 
  select (Location, Dia, h_ini, individual, plot) %>% 
  group_by (Location, Dia, plot, individual) %>% 
  summarise (h_ini = first (h_ini))
## agafo la relació d'individus detectats en cada dia, lloc, plot i hora de mostreig


sum.fem <- read_excel ("C:/Users/MariaCyborg/OneDrive - CREAF/Dades i analisis/Bases definitives/Behaviour_oviposition/Behaviour_20172018_20190723.xlsx",
                       sheet = 5) %>% 
  mutate (Date = as_date (Date),
          Code = as.character (Code),
          month = month (Date)) %>%
  filter (month >= 6,
          month <= 9) %>% 
  left_join (fems, by = c ("Location" = "Location",
                           "Date" = "Dia",
                           "Code" = "individual")) %>% 
  select (-Land_num, -COMMENTS, )
## afegeixo el plot i l'hora de mostreig a la relació de femelles detectades i el nombre d'ous que han post

fem.pois1 <- read_excel ("C:/Users/MariaCyborg/OneDrive - CREAF/Dades i analisis/Bases definitives/Behaviour_oviposition/Behaviour_20172018_20190723.xlsx",
                         sheet = 4) %>% 
  mutate_if (is.character, as.factor) %>% 
  rename (duration = "Duration (min)",
          rate_ovi = "Rate_ovi (obs/h)",
          rate_ind_ovi = "Rate_ind_ovi (ind/h)",
          outlier = "Outlier (1=yes)") %>% 
  mutate (initial_h = as_hms (Intial_H),
          final_h = as_hms (Final_H),
          day_period = case_when (Time_zone < 11 ~ "1",
                                  Time_zone < 13 ~ "2",
                                  Time_zone < 15 ~ "3",
                                  Time_zone < 17 ~ "4",
                                  Time_zone >= 17 ~ "5"),
          day_period = as.factor (day_period),
          Habitat_class = factor (Habitat_class, levels = c ("C", "OC", "O")),
          Date = as_date (Date)) %>%  
  filter (Month >= 6,
          Month <= 9) %>% 
  left_join (sum.fem, by = c("Location" = "Location",
                             "Date" = "Date",
                             "Habitat_class" = "Land",
                             "Species" = "Species",
                             "initial_h" = "h_ini",
                             "Plot" = "plot",
                             "Month" = "month")) %>%
  group_by (Location, Date, Year, Month, Habitat_class, Species, day_period, Code) %>% 
  summarise (n_obs_real = sum (Numb_Ovi, na.rm = T),
             n_fem_real = sum (Ind_ovi, na.rm = T),
             n_obs_estm = sum (rate_ovi, na.rm = T)*2,
             n_fem_estm = sum (rate_ind_ovi, na.rm = T)*2,
             duration = sum (duration, na.rm = T),
             n_plots = n_distinct (Plot),
             Num_obs = sum (Num_obs, na.rm = T)) %>% 
  mutate (Code = case_when ((n_fem_real > 0 & is.na (Code) & Species == "PN") ~ "627",
                            (n_fem_real > 0 & is.na (Code) & Species == "PR") ~ "629",
                            TRUE ~ as.character (Code)),
          Num_obs = case_when (Code == "627" ~ 2,
                               Code == "629" ~ 1,
                               is.na (Code) ~ 0,
                               TRUE ~ as.numeric (Num_obs)))

## a cada mostreig fet, li afegeixo les diverses femelles detectades. Ara, cada entrada és una femella diferent.
## Problemes:
# femelles que no tenien codi (627, 628 i 629)
# en la meva separació per mostrejos (surveys), no distingeixo plots sinó microambients, però en fer el merging, tot i ja no tenir el camp plot, el nombre de femelles i d'observacions encara està per plot, i no per microambient.
# si en algun plot d'algun microambient no s'havien vist femelles, en el summary per microambient se'm manté una fila de codi absent.

missing.case <- fem.pois1 %>% 
  filter (Code == "627") %>% 
  mutate (Code = "628",
          Num_obs = 1)

fem.pois2 <- fem.pois1 %>% 
  full_join (missing.case) %>% 
  group_by (Location, Date, Year, Month, Habitat_class, day_period) %>% 
  mutate (survey = group_indices()) %>% #afegeixo survey
  group_by (survey, Species) %>% 
  mutate (check.na = n_distinct (Code))


fem.pois <- fem.pois2 %>% 
  distinct (survey, Species, n_fem_real, n_obs_real, duration) %>% 
  group_by (survey, Species) %>% 
  summarise (n_fem_real = sum (n_fem_real, na.rm = T),
             n_obs_real = sum (n_obs_real, na.rm = T),
             duration = sum (duration, na.rm = T)) %>% 
  right_join (select(fem.pois2, -c(n_fem_real, n_obs_real, duration))) %>% 
  filter ((!is.na (Code) & n_fem_real > 0) |
            is.na (Code) & n_fem_real == 0) %>% 
  rename (n_fem_surv = n_fem_real,
          n_obs_surv = n_obs_real) %>% 
  select (survey, Location, Date, Year, Month, day_period, Habitat_class, Species,
          everything (),
          -c(n_fem_estm, n_obs_estm, check.na)) %>% 
  rename (site = Location,
          date = Date,
          year = Year,
          month = Month,
          microhabitat = Habitat_class,
          species = Species,
          fem_code = Code,
          n_ovi = Num_obs) %>% 
  mutate (jday = yday (date),
          jday_cat = as.character (jday)) %>% 
  select (survey, site, date, year, month, jday, jday_cat, day_period:species, duration, fem_code, n_ovi)

rm (missing.case, sum.fem, fem.pois1, fem.pois2, fems)


# n_fem/survey (model 5.4, posat al supp)
#########################################

#no permet distingir els ous de cada femella

surv <- fem.pois %>% 
  select (survey, site, jday, day_period, microhabitat) %>% 
  distinct ()


ovi.pois <- read_excel ("C:/Users/MariaCyborg/OneDrive - CREAF/Dades i analisis/Bases definitives/Behaviour_oviposition/Behaviour_20172018_20190723.xlsx",
                        sheet = 4) %>% 
  mutate_if (is.character, as.factor) %>% 
  rename (duration = "Duration (min)",
          rate_ovi = "Rate_ovi (obs/h)",
          rate_ind_ovi = "Rate_ind_ovi (ind/h)",
          outlier = "Outlier (1=yes)") %>% 
  mutate (initial_h = as.hms (Intial_H),
          final_h = as.hms (Final_H),
          day_period = case_when (Time_zone < 11 ~ "1",
                                  Time_zone < 13 ~ "2",
                                  Time_zone < 15 ~ "3",
                                  Time_zone < 17 ~ "4",
                                  Time_zone >= 17 ~ "5"),
          day_period = as.factor (day_period),
          Habitat_class = factor (Habitat_class, levels = c ("C", "OC", "O"))) %>% 
  group_by (Location, Date, Year, Month, Habitat_class, Species, day_period) %>% 
  summarise (n_obs_real = sum (Numb_Ovi, na.rm = T),
             n_fem_real = sum (Ind_ovi, na.rm = T),
             n_obs_estm = sum (rate_ovi, na.rm = T)*2,
             n_fem_estm = sum (rate_ind_ovi, na.rm = T)*2,
             duration = sum (duration, na.rm = T),
             n_plots = n_distinct (Plot)) %>% 
  rownames_to_column (var = "survey") %>% 
  rename (site = Location,
          date = Date,
          year = Year,
          month = Month,
          microhabitat = Habitat_class,
          species = Species,
          n_ovi = n_obs_real,
          n_fem = n_fem_real) %>% 
  mutate (jday = yday (date),
          jday_cat = as.character (jday)) %>% 
  filter (month >= 6, month <= 9) %>% 
  select (site, date, year, month, jday, jday_cat, microhabitat:n_fem, duration) %>% 
  left_join (surv) %>% 
  select (survey, everything())


# temperatures oviposicio
#########################

micro_period <- fem.pois %>% 
  select (fem_code, day_period, microhabitat) %>% 
  filter (!is.na (fem_code))

ovi.temp <- read_excel ("C:/Users/MariaCyborg/OneDrive - CREAF/Dades i analisis/Bases definitives/Behaviour_oviposition/Behaviour_20172018_20190723.xlsx",
                        sheet = 2) %>% 
  mutate (Dia = as.numeric (Dia),
          Dia = as_date (Dia, origin = "1899-12-30")) %>% 
  rename (h_ini = "H inicial",
          h_fin = "H final",
          species = "Espècie papallona",
          fem_code = "Número individu",
          plot = "Parcel·la",
          T_upp = "Tª anvers",
          T_und = "Tª revers") %>% 
  mutate (h_ini = as_hms (h_ini),
          h_fin = as.numeric (h_fin),
          h_fin = h_fin*60*60*24,
          h_fin = as_hms (h_fin),
          jday = yday (Dia),
          jday_cat = as.character (jday)) %>% 
  filter (COMP == "O",
          (species == "PN" | species == "PR"),
          (Month == 6 | Month == 7 | Month == 8 | Month == 9)) %>% 
  mutate (fem_code = case_when ((fem_code == 376 & T_upp > 27) ~ "627",
                                fem_code == 376 ~ "628",
                                fem_code == 507 ~ "629",
                                TRUE ~ fem_code)) %>% 
  select (Location, Dia, YEAR, Month, jday, jday_cat, species, fem_code, T_upp, T_und) %>% 
  left_join (micro_period) %>% 
  select (survey, Location, Dia, YEAR, Month, jday, jday_cat, day_period, microhabitat, species, fem_code, T_upp, T_und) %>% 
  rename (site = Location,
          date = Dia,
          year = YEAR,
          month = Month) %>% 
  filter_at (vars (T_upp, T_und), any_vars (!is.na (.)))

rm(surv, micro_period)

write_csv (fem.pois, "final_versions/ecotones_procB/ovi_by_female.csv")
write_csv (ovi.pois, "final_versions/ecotones_procB/ovi_by_survey.csv")
write_csv (ovi.temp, "final_versions/ecotones_procB/ovi_temp.csv")

########################
# Microclimatic regime #
########################

sensor <- read.csv ("C:/Users/MariaCyborg/OneDrive - CREAF/Dades i analisis/Bases definitives/Thermal_environment/SENSORS_AE_CJ_2018_v4.csv", dec = ",", sep = ";")

vari.sens <- sensor %>% 
  mutate (patch_plot =
            case_when ((patch == "LD41" |
                          patch == "BN3" |
                          Sensor_name == "s12_alzinar_pal_telefon") ~ "O",
                       (patch == "CJ4" |
                          patch == "CJ3_2" |
                          patch == "BN41" |
                          patch == "BN2") ~ "OC",
                       (patch == "LD2" |
                          patch == "LD1" |
                          patch == "BN42" |
                          patch == "CJ2") ~ "C",
                       TRUE ~ NA_character_),
          patch_plot = factor (patch_plot, levels = c ("C", "OC", "O"))) %>% 
  filter (!is.na (patch_plot),
          (month == 6 | month == 7 | month == 8),
          Sensor_name != "pont_brassica_s6") %>% #on estava exactament el sensor?
  mutate (Sensor_name = as.character (Sensor_name),
          Sensor_name = as.factor (Sensor_name),
          winter_date = as.character (winter_date),
          winter_date = parse_date_time(winter_date, "%d/%m/%Y", exact = T)) %>% 
  select (site, patch_plot, Sensor_name, winter_date, winter_jday:winter_year, TEMP:DP) %>% 
  rename (microhabitat = patch_plot)

metcat <- read_excel ("C:/Users/MariaCyborg/OneDrive - CREAF/Dades i analisis/Bases definitives/Thermal_environment/MeteoCat_AE_CJ.xlsx",
                      col_types = c(rep ("guess", 13),
                                    "numeric", "numeric",
                                    "skip", "skip", "skip")) %>% 
  mutate_if (is.character, as.factor) %>% 
  filter (Year <= 2017, Year >= 2014) %>% 
  rename (site = Site,
          year = Year,
          jday = "Julian day",
          METCAT_TM = TM,
          METCAT_TX = TX,
          METCAT_TN = TN) %>% 
  mutate (site = recode (site, "AIGUAMOLLS EMPORDA" = "AE", "CAN JORDA" = "CJ"))

vari.sens.sum <- vari.sens %>% 
  group_by (site, microhabitat, Sensor_name, winter_date) %>% 
  mutate (count = n()) %>% 
  filter (count == 24) %>% 
  summarise (mean_t = mean (TEMP, na.rm = T),
             min_t = min (TEMP, na.rm = T),
             max_t = max (TEMP, na.rm = T),
             range = max_t - min_t,
             sd = sd (TEMP, na.rm = T),
             var = var (TEMP, na.rm = T),
             cv = sd/mean_t) %>% 
  left_join (metcat,
             by = c("winter_date" = "Date", "site" = "site")) %>% 
  mutate (mean_amp = mean_t - METCAT_TM,
          max_amp = max_t - METCAT_TX,
          min_amp = min_t - METCAT_TN,
          species = case_when (microhabitat == "OC" ~ "P. napi",
                               microhabitat == "O" ~ "P. rapae",
                               TRUE ~ NA_character_)) %>% 
  select (site:winter_date, year, Season, Month, jday, mean_t:cv, METCAT_TM:HRM, contains("amp"), species)

vari.sens.sum.mc <- metcat %>% 
  rename (mean_t = METCAT_TM,
          max_t = METCAT_TX,
          min_t = METCAT_TN,
          winter_date = Date,
          hrm = HRM) %>% 
  mutate (Sensor_name = "Meteo",
          microhabitat = "M",
          range = max_t-min_t,
          min_hr = NA) %>% 
  select (site, winter_date, Sensor_name, microhabitat,
          mean_t, max_t, min_t, range, hrm, min_hr) %>% 
  right_join (vari.sens.sum, by = c("site" = "site",
                                    "winter_date" = "winter_date"),
              suffix = c ("", "_S")) %>% 
  select (-ends_with ("_S"), -(sd:cv)) %>% 
  distinct () %>% 
  bind_rows (vari.sens.sum) %>% 
  mutate (microhabitat = factor (microhabitat,
                               levels = c ("C", "OC", "O", "M")),
          microhabitat_2 = factor (microhabitat,
                                 levels = c ("C", "OC", "O"))) %>% 
  select (site, microhabitat, Sensor_name, winter_date, year:jday, mean_t:min_hr, species)

write_csv (vari.sens, "final_versions/ecotones_procB/sensor_hourly.csv")
write_csv (vari.sens.sum, "final_versions/ecotones_procB/sensor_daily.csv")
write_csv (vari.sens.sum.mc, "final_versions/ecotones_procB/sensor_daily_mc_as_micro.csv")


########
# TDT? #
########