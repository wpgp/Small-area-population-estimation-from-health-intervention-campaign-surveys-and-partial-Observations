####--TITLE:Small area population estimation from health intervention campaign surveys and partially observed settlement data--#####
####--METHODS: GEOSTATISTICAL BAYESIAN HIERARCHICAL REGRESSION MODEL
####--AUTHOR: DR CHIBUZOR CHRISTOPHER NNANATU
####--INSTITUTION: WORLDPOP, UNIVERSITY OF SOUTHAMPTON 
####--DATE:19 DECEMBER 2022. REVIEWED: 25/03/2024
###=====================================================================================================
#           
#      DATA PREPARATION AND COVARIATES SELECTION
#      -----------------------------------------
#------------------------------------------------------------------------------------------------

rm(list=ls()) #---Clear workspace

##---List of key packages
packages <- c("raster", "haven", "sf","sp", "tidyverse","rgdal",
              "lattice", "gridExtra", "devtools", "rlang", "spdep", "viridis", "MASS")

##----Only install packages listed and not yet installed
if(length(setdiff(packages, rownames(installed.packages()))) > 0) { 
  install.packages(setdiff(packages, rownames(installed.packages())), type="binary") }

#Installing INLA!!
if(length(setdiff("INLA", rownames(installed.packages()))) > 0){
  install.packages("INLA", type="binary", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
}

#install.packages("INLA")
library(INLA)
lapply(packages, library, character.only = TRUE) ##--access the libraries

##-------------------------------------------------------------------------------
###--LOAD DATA FROM GITHUB REPOSITORY
#--------------------------------------------------------------------------------
# Demographic data: CU - level survey data
surv_df <- read.csv("https://raw.github.com/wpgp/Small-area-population-estimation-from-health-intervention-campaign-surveys-and-partial-Observations/main/surv_data.csv")
names(surv_df)
str(surv_df)




#---------------------------------------------------------------------------
# covariates Z-score scaling 
stdize <- function(x)
{
  stdz <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  return(stdz)
}

#--------------------------------------------------------------------------
# More data preparation s
#----------------------------------------------------------------------------
datt <- surv_df #-------duplicate data
names(datt)
mod_covs <- 10:61 ; length(10:61)#----52 covariates


head(datt[, mod_covs]) #-----view the first sic rows of the covariates columns

# recode covariates (this is just for convenience)
colnames(datt)[mod_covs] <- paste0("x", 1:52, sep="") #--rename covariates 
names(datt)  #------check

# apply covariates scaling to the model covariates only
datt[,mod_covs] <- apply(datt[,mod_covs], 2, stdize)   


#-----------------------------------------------------------------------
#  POPN: Number of people per CENSUS UNIT (CU)
#  BLDG21: aggregated building intensity per CU
#  TYPE: Settlement type
#-----------------------------------------------------------------------

# explore the datasets 
par(mfrow=c(2,2))
hist(datt$POPN)  #---observed pop total per CU (right skewed)
boxplot(datt$POPN) # -- one cu with more than 8000 people
hist(datt$BLDG21) #---Building intensity  (right skewed)


par(mfrow=c(1,1)) # restore to one image window


# Subset data into training and test sets
dim(dat_train <- datt[!is.na(datt$POPN),]) #-- Training set - all CUs with Observations only (16872 CUs)
dim(dat_pred <- datt[is.na(datt$POPN),])   #-- Prediction set - all unsampled CUs (15095 CUs)

#----------------Covariates Selection-----------------------------------
# Run GLM-based stepwise selection 
# a few more relevant libraries
library(car) ##--For calculating variance inflation factor (vif)
library(dplyr)
library(tidyverse)


names(dat_train)
dat_train$resp <- dat_train$POPN # create the 'resp' variable

covs_pop <- dat_train[,c(62,mod_covs)] #--response of interest is found in column 62
covs_pop # check


covs_pop <- covs_pop %>% drop_na() #-- drop all NAs
dim(covs_pop <- covs_pop[is.finite(covs_pop$resp),]) #--remove all infinite values (if any)

#- fit a GLM based model 
fpop <- glm.nb(resp ~., data=covs_pop) #---negative binomial
                                      # note: continous population density will use either Gaussian or Gamma density
summary(fpop)


# run stepwise regression using 'both' forward and backward algorithms
step_pop <- stepAIC(fpop, scale = 0,
                    direction = c("both"),
                    trace = 1, keep = NULL, steps = 1000, use.start = FALSE,
                    k = 2)
step_pop


# select the best model
fpop2 <- glm.nb(formula = resp ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + 
         x9 + x10 + x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + 
         x19 + x20 + x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + 
         x29 + x30 + x31 + x32 + x33 + x34 + x35 + x36 + x37 + x38 + 
         x39 + x40 + x41 + x42 + x43 + x45 + x46 + x47 + x48 + x49 + 
         x50 + x51 + x52, data = covs_pop, init.theta = 1.889859884, 
       link = log)

summary(fpop2)

# only sttaistocally significant variables with vif less than 5 are retained
vif_pop = vif(fpop2)
vif_pop[which(vif_pop < 5)]


# Refit model with only covariates with vif values less than 5
fpop3 <- glm.nb(formula = resp ~x3 +  x11 + x15 + x29 + x30 + x31 + x32 + x33  + x35 + x36  + x38 + 
                  x39 +  x45 + x46 + x48 + x51 + x52, data = covs_pop, init.theta = 1.889859884, 
                link = log)

summary(fpop3)


# Drop variables not statistically significant and refit
fpop4 <- glm.nb(formula = resp ~x3 +  x11 + x15 + x29 + x30 + x31 + x32 + x33  + x35 + x36 + # only 'x38' was dropped
                  x39 +  x45 + x46 + x48 + x51 + x52, data = covs_pop, init.theta = 1.889859884, 
                link = log)

summary(fpop4) 


# Drop variables not statistically significant and refit
fpop5 <- glm.nb(formula = resp ~x3 +  x11 + x15 + x29 + x30 + x31 + x32 + x35 + x36 + # only 'x33' was dropped
                  x39 +  x45 + x46 + x48 + x51 + x52, data = covs_pop, init.theta = 1.889859884, 
                link = log)

summary(fpop5) # all now statistically significant with p-value less than 0.05



# recalculate vif to ensure effects of multicollinearity is drastically reduced
vif_pop5 = vif(fpop5)
vif_pop5[which(vif_pop5 < 5)]



# Best fit covariates eventually selected
cov_pop <- c("x3","x11","x15","x29","x30", "x31","x32","x35", "x36",
             "x39","x45","x46","x48","x51","x52")

cov_names <- names(dat_train)[mod_covs]
covariates_pop <- cov_names[c(3, 11, 15, 29, 30, 31, 32, 
                              35, 36, 39, 45, 46, 48, 51, 52)]
pop_covs <- data.frame(cov = cov_pop, name=covariates_pop)
#-------------------------------------

##----Make correlation plot
#install.packages("corrplot")
require(corrplot)
#png(paste0(results_path,"cor_plots.png"))
corrplot(
  cor(covs_pop[,cov_pop]),
  method = 'square',
  type = 'upper',
  tl.col = 'black',
  tl.cex = 1,
  col = colorRampPalette(c('purple', 'dark green'))(200)
)
#dev.off()


