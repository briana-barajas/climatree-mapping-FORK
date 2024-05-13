# Mapping Repository
Repository for data wrangling and spatial analysis

Definitions of key terms:

**sensitivity:** The partial effect of CWD on tree growth (the direction and degree to which trees will react to changes in water availability)

**vulnerability:** Combining sensitivity values with predictions of CWD change (using CMIP5) to predict how climate change will effect future tree growth. 

**Coupled Model Intercomparison Project Phase 5 (CMIP5):** An ensemble of climate prediction models used to predict CWD change between historic (1970-2000) and end-of-century (2090-2100)  

## This repository has 5 Rmarkdown files in 2 categories:

###[Category 1:] First 2 files that were modified to run with updated TerraClimate data. The outputs from these files were used as inputs in the Category 2 files below.

**1_climate_niche.R** --> This markdown document explores how historical climate impacts weather sensitivity of tree species. 

**2_plot_level_regression.R** --> This document uses plot-level regressions of RWI sensitivity to annual weather variations. 

###[Category 2:] 3 files of scripts modified from a global to a species-level analysis of tree sensitivity and vulnerability.

**3_run_regressions.R** --> This script species-level regressions to analyze the impact of historical climate on weather sensitivity (combines CMIP5 data with TerraClimate data)

**4_sens_predictions.R** --> This script creates predictions of how climate change will impact growth across each species range (vulnerability). 

**5_mapping.R** --> This script creates sensitivity and vulnerability maps for targeted species. 


