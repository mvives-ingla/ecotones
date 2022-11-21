Data and code associated with the manuscript  
Vives-Ingla, M., J. Sala-Garcia, C. Stefanescu, A. Casadó-Tortosa, M. Garcia, J. Peñuelas, and J. Carnicer. 2022. Interspecific differences in microhabitat use expose insects to contrasting thermal mortality. Ecological Monographs.

<b>contact:</b>  
m.vives[at]creaf.uab.cat  
jofre.carnicer[at]ub.edu


---
## Abstract

Ecotones linking open and forested habitats contain multiple microhabitats with varying vegetal structure and microclimatic regimes. Ecotones host many insect species whose development is intimately linked to the microclimatic conditions where they grow (e.g. the leaves of their host plants and the surrounding air). Yet microclimatic heterogeneity at these fine scales and its effects on insects remain poorly quantified for most species. Here we studied how interspecific differences in the use of microhabitats across ecotones lead to contrasting thermal exposure and survival costs between two closely-related butterflies (<i>Pieris napi</i> and <i>P. rapae</i>). We first assessed whether butterflies selected different microhabitats to oviposit and quantified the thermal conditions at the microhabitat and foliar scales. We also assessed concurrent changes in the quality and availability of host plants. Finally, we quantified larval time of death under different experimental temperatures (TDT curves) to predict their thermal mortality considering both the intensity and the duration of the microclimatic heat challenges in the field. We identified six processes determining larval thermal exposure at fine scales associated with butterfly oviposition behavior, canopy shading and heat and water fluxes at the soil and foliar levels. Leaves in open microhabitats could reach temperatures 3–10 ºC warmer than the surrounding air while more closed microhabitats presented more buffered and homogeneous temperatures. Interspecific differences in microhabitat use matched the TDT curves and the thermal mortality in the field. Open microhabitats posed acute heat challenges that were better withstood by the thermotolerant butterfly, <i>P. rapae</i>, where the species mainly laid their eggs. Despite being more thermosensitive, <i>P. napi</i> was predicted to present higher survivals than <i>P. rapae</i> due to the thermal buffering provided by their selected microhabitats. However, its offspring could be more vulnerable to host-plant scarcity during summer drought periods. Overall, the different interaction of the butterflies with microclimatic and host plant variation emerging at fine scales and their different thermal sensitivity posed them contrasting heat and resource challenges. Our results contribute to set a new framework that predicts insect vulnerability to climate change based on their thermal sensitivity and the intensity, duration and accumulation of heat exposure.


---
## __`data/`__
 
  - __`but_pheno.csv`__: weekly counts of butterflies in the study sites during 2017    
  - __`data_fruits.csv`__: number of fruits per host-plant species and microhabitat  
  - __`data_leaves.csv`__, __`data_microhabitat.csv`__, __`data_plants.csv`__: microclimatic and host-plant traits monitored in 2017 at the foliar, microhabitat and plant scales.  
  - __`sensors_hourly.csv`__: hourly records with LASCAR data loggers at the microhabitat scale
  - __`thermal_gradient_height.csv`__: thermal profiles from soil surface to 1m height
  - __`meteorological_data.csv`__: data from nearby standardised weather stations  
  - __`ovi_*`__: oviposition data (an entry per ovipositing female, an entry per census, and temperatures recorded during oviposition)  
  - __`TDT_experiment.txt`__: thermal death times of the larvae in heat tolerance essays


---
## __`code/`__
  - __`figs/`__: code that reproduce all the figures contained in the manuscript  
  - __`thermal_mortality/`__:   
    - __`sim_thermal_mortality.R`__: code to predict thermal mortality from microclimatic data and bootstrapped TDT curves.  
    - __`Thermal_landscape_functions_mod.R`__: code with the functions needed to predict thermal mortality in dynamic thermal conditions from a TDT curve. Modified from Rezende, E., F. Bozinovic, A. Szilágyi, and M. Santos. 2020. Dataset and scripts from: Predicting temperature mortality and selection in natural Drosophila populations.  
    - __`function_trunk_mortality`__: custom functions to predict thermal mortality from a tolerance landscape and a temperature series for an organism that activates thermal avoidance behaviors when a thermal threshold is exceeded (simulated by trunkating the thermal profile at this thershold).  
  - __`*_models`__: code with the (G)LMMs included in the analyses.

  
