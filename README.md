Data and code associated with the manuscript  
Vives-Ingla, M., J. Sala-Garcia, C. Stefanescu, A. Casadó-Tortosa, M. Garcia, J. Peñuelas, and J. Carnicer.  
Interspecific differences in microhabitat use expose insects to contrasting thermal mortality.

+ code/  
  - figs/: code that reproduce all the figures contained in the manuscript  
  - thermal_mortality/:   
    - sim_thermal_mortality.R: code to predict thermal mortality from microclimatic data and bootstrapped TDT curves.  
    - Thermal_landscape_functions_mod.R: code with the functions needed to predict thermal mortality in dynamic thermal conditions from a TDT curve. Modified from Rezende, E., F. Bozinovic, A. Szilágyi, and M. Santos. 2020. Dataset and scripts from: Predicting temperature mortality and selection in natural Drosophila populations.  
    - function_trunk_mortality: functions to predict thermal mortality from a tolerance landscape and a thermal profile for an organism that activates thermal avoidance behaviors when a thermal threshold is exceeded (simulated by trunkating the thermal profile at this thershold).  
  - *_models: code for the (G)LMMs included in the analyses.
+ data/  
  - but_pheno.csv: weekly counts of butterflies in the study sites    
  - data_fruits.csv: number of fruits per species and microhabitat  
  - data_leaves.csv, data_microhabitat.csv, data_plants.csv: microclimatic and host-plant traits monitored in 2017 at the foliar, microhabitat and plant scales.  
  - meteorological_data.csv: data from weather stations  
  - ovi_*: oviposition data (an entry per ovipositing female, an entry per census, and temperatures recorded during oviposition)  
  - TDT_experiment.txt: thermal death times of the larvae in heat tolerance essays
  
