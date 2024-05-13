# Mapping Repository 🌲
Repository for data wrangling and spatial analysis

* to access original scripts before modification by this project, vist Treeconomics GitHub Repository: https://github.com/rheilmayr/Treeconomics 

## Data and Term Definitions

<u>Data:<u>

Raw data file "essential_cwd.csv" can be accessed from the International Tree Ring Data Bank (ITRDB): https://www.ncei.noaa.gov/products/paleoclimatology/tree-ring

Raw climate data files "TerraClimate19611990_pet.nc" and "TerraClimate19611990_def.nc" can be accessed from the TerraClimate website: https://www.climatologylab.org/terraclimate.html

* For access to other data sources, contact Dr. Joan Dudney 

Definitions of key terms:

**CWD(climatic water deficit):** A measure of how much less water availability there will be in the future

**PET(potential evapotranspiration):** The maximum amount of growth a tree will have, when water is unlimited

**Sensitivity:** The partial effect of CWD on tree growth (the direction and degree to which trees will react to changes in water availability)

**Vulnerability:** Combining sensitivity values with predictions of CWD change (using CMIP5) to predict how climate change will effect future tree growth. 

**Coupled Model Intercomparison Project Phase 5 (CMIP5):** An ensemble of climate prediction models used to predict CWD change between historic (1970-2000) and end-of-century (2090-2100)  

## This repository has 5 Rmarkdown files:

***1_climate_niche.R** --> This markdown document explores how historical climate impacts weather sensitivity of tree species. 

**2_plot_level_regression.R** --> This document uses plot-level regressions of RWI sensitivity to annual weather variations. 

**3_run_regressions.R** --> This script species-level regressions to analyze the impact of historical climate on weather sensitivity (combines CMIP5 data with TerraClimate data)

**4_sens_predictions.R** --> This script creates predictions of how climate change will impact growth across each species range (vulnerability). 

**5_mapping.R** --> This script creates sensitivity and vulnerability maps for targeted species. 

## Repository Structure:
```
climatree-mapping-repo
├── main.R
├── 1_climate_niche.R
├── 2_plot_level_regressions.R
├── 3_run_regressions.R
├── 4_sens_predictions.R
├── 5_mapping.R
└── prep_scripts
    ├── 2_spp_list.R
    └── create_top_species.R

```
