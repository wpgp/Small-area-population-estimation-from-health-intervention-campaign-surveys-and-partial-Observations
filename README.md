# Small-area-population-estimation-from-health-intervention-campaign-surveys-and-partial-Observations - Nnanatu et al. (2024)

## Background
This repository provides the R scripts for the implementation of the statistical modelling appraoch described within the 'Small-area-population-estimation-from-health-intervention-campaign-surveys-and-partial-Observations' paper. It also includes the key (aggregated and anonymized) input datasets used for the modelling. Please note that the datasets are included here for learning purposes only. Additionally, the simulation study R scripts are also included. 
This ReadMe file contains names and the descriptions of the various R codes and the associated datasets to facilitate learning and reproducibility. The outlines are provided for ease of use by both beginners and advanced R programming language users. 

Below, we provide the details of system requirements and software installation guides followed by the the descriptions of the R files and the associated data files. 

## System Requirements
The R codes provided within this repository can be run within the operating system (e.g., Windows, Mac) of any machine with R software downloaded and installed. You can either use the latest R version or if you need to reproduce the model estimates, please kindly use R version 4.0.2. There would be slight differences in the parameter estimates when different versions of R are used.

## Installation Guide for R, RStudio and R packages
### To install R on your machine (e.g., Windows Operating System), please 
1) visit CRAN website (https://cran.r-project.org/)
2) click on "Download R for Windows"
3) click on "Install R for the first time": This allows you to download the R excutable file or .exe file
4) Run the R executable file and follow the on-screen instructions on your screen and allow the app to make changes to your device.
5) Choose the installation language and follow the on-screen instructions on your screen
6) When completed, click on 'finish' and you are ready to launch R

### To install RStudio on your machine (e.g., Windows Operating System), please 
To download and install RStudio, do the following:
1) Go to https://posit.co/download/rstudio-desktop/
2) Click on "Download RStudio Desktop"
3) Select the recommended version and save the executable file
4) Run executable file while following the on-screen guide

### To install R packages on your machine after the installation of R
use the command 'install.package("name of the package")' followed by 'library("name of the package")' to have access to the libarary and depedencies. 

## Demo
A  brief demostration of how to install an R package is provided below: 
Assuming that I would like ton install the tidyverse package, I would type the following commands in R script window:
install.packages("tidyverse")

library(tidyverse)

## Instructions for use
To use the scripts provided within this repository, you need to make sure you 
1) Already have R installed and open
2) On the left pane of the repository, locate and click on the R file of interest
3) As soon as the file opens on your GitHub, you can copy the codes and paste onto the R script on your computer
4) Although all the input datasets are mapped onto the GitHub repository which means that you do not have to download and save them locally, you may require to create an output folder to assess the results locally on your machine.

## Scripts and input data files descriptions 
### R scripts
There are four R scripts included within this folder: 
1) covariates_selection.R - this was to select the best fit covariates that woulc be used for the final model. A stepwise regression variable selection techniques with forward and backward (or 'both')algrorithms.
2) small_area_pop_estimation_simulation_study.R - this script covers the simulation stusy we conducted to assess the sensitivity of the BHM (traditional) and TSBHM (proposed) approaches over different parameter values/levels of missingness in the input datasets.
3) small_area_pop_estimation_sim_study_graphs.R - This is used to produce the graphs cited in the main document under simulation study.
4) small_area_pop_estimation_application.R - this script outlines the implementation of our methodology using real data collected across the various census units in PNG

### Input Data
Below are two major input datasets used: 
1) survey_data.RData: this is the main input data containing counts of people, building counts, and all the geospatial covariates prepared at the census unit level. (please use only 'survey_data.RData'. It contains a data frame called 'covs' which ais the main data. 
2) cu_boundary.gpkg: this matches perfectly with the surv_data and contains most of the variables in the .csv file. The centroids of the c_boundary file is used as the longitude and latituide for the surv_data file. 

## Other files
We have included the supplementary document containing extra information that could not go into the main manuscript here. The file is called "Supplemental-Small-Area_Pop_Estimates.pdf"

#### Please kindly leave us a feedback to help improve our future outputs. You can contact the corresponding author on cc.nnanatu@soton.ac.uk or nnanatuchibuzor@gmail.com - thanks! 
