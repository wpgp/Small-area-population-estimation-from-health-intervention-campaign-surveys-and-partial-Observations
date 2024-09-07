# Small-area-population-estimation-from-health-intervention-campaign-surveys-and-partial-Observations
This repository provides the R scripts for the 'Small area population estimation from health intervention data ...' paper. The scripts include the simulation study and application studies.  

Below, the description of each of the files included within thos repository is provided:

R scripts
There are four R scripts included within this folder - these are 
1) covariates_selection.R - this was to select the best fit covariates that woulc be used for the final model. A stepwise regression variable selection techniques with forward and backward (or 'both')algrorithms.
2) small_area_pop_estimation_simulation_study.R - this script covers the simulation stusy we conducted to assess the sensitivity of the BHM (traditional) and TSBHM (proposed) approaches over different parameter values/levels of missingness in the input datasets.
3) small_area_pop_estimation_sim_study_graphs.R - This is used to produce the graphs cited in the main document under simulation study.
4) small_area_pop_estimation_application.R - this script outlines the implementation of our methodology using real data collected across the various census units in PNG

Input Data
1) surv_data: this is the input data containing counts of people, building counts, and all the geospatial covariates prepared at the census unit level. (there are two versions of the file - text or .txt (surv_data.txt) and .csv (surv_data.csv)
2) cu_boundary.gpkg: this matches perfectly with the surv_data and contains most of the variables in the .csv file. The centroids of the c_boundary file is used as the longitude and latituide for the surv_data file. 
