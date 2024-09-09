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
#Demographic data: Census unit (CU) - level survey data (contains observed count of people per CU)
githublink <- "https://raw.github.com/wpgp/Small-area-population-estimation-from-health-intervention-campaign-surveys-and-partial-Observations/main/survey_data.RData"
load(url(githublink))
names(covs) # the data frame is called 'covs' 

# CU level shapefile
cu_shp <- st_read("https://raw.github.com/wpgp/Small-area-population-estimation-from-health-intervention-campaign-surveys-and-partial-Observations/main/cu_boundary.gpkg")
names(cu_shp)
str(cu_shp)
plot(cu_shp["Shape_Area"])


# convert sf shapefile to sp 
shp <- as(st_geometry(cu_shp), "Spatial") # converts sf data to sp. 
plot(shp)

##
rm(cu_shp) # remove duplicated filr to conserve memory

###
ls() # view console contents
covDF <- as.data.frame(covs) # convert to a data frame


#plot(covDF$lon, covDF$lat)
#names(covDF); str(covDF); dim(covDF); dim(shp)# shp has one extra row

shp <- shp[-32100,] #--remove row 32100 from the shapefile as it does not exist
dim(covDF); dim(shp) #---confirm
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
##----Extract the coordiates - the centroids of the CUs
#--------------------------------------------------------------------
library(sp)
covDF$lon <- coordinates(shp)[,1]#--add lon-lat to the data
covDF$lat <- coordinates(shp)[,2]
datt <- covDF
#---------------------------------------------------------------------------



#---------------------------------------------------------------------------
# covariates Z-score scaling function
stdize <- function(x)
{
  stdz <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  return(stdz)
}

#-----------------------------------------------------------------------------
# Function: Model fit metrics 
mod_metrics2 <- function(obs, pred)
{
  residual = pred - obs
  INACCURACY = mean(abs(residual), na.rm=T)#MAE
  MSE = mean(residual^2, na.rm=T)
  RMSE = sqrt(MSE)
  BIAS = mean(residual, na.rm=T)
  IMPRECISION = sd(residual, na.rm=T)
  COR = cor(obs[!is.na(obs)], pred[!is.na(obs)])
  output <- list(MAE  = INACCURACY ,
                 RMSE = RMSE,
                 BIAS = abs(BIAS),
                 COR=COR)
  return(output)
}

#--------------------------------------------------------------------------
# More data preparation s
#----------------------------------------------------------------------------


#---------------------------------------------------------------------
#           BAYESIAN STATISTICAL MODELLING
#-------------------------------------------------------------------
# Define rhe points for both training and testint/predictions separately
# Subset data into training and test sets
dim(dat_train <- datt[!is.na(datt$POPN),]) #-- Training set - all CUs with Observations only (16872 CUs)
dim(dat_pred <- datt[is.na(datt$POPN),])   #-- Prediction set - all unsampled CUs (15095 CUs)


# specify GPS coordinates for each dataset  
coords = cbind(datt$lon, datt$lat)
coords_train = cbind(dat_train$lon, dat_train$lat)
coords_pred = cbind(dat_pred$lon, dat_pred$lat)

##-View the study locations
plot(coords, cex=0.1, col="black")
plot(coords_train, cex=0.1, col="blue",
     ylab="Latitude",
     xlab="Longitude",
     cex.axis=1.5)
points(coords_pred, cex = 0.1, col="red")


###----------------------------
# Further explorations of the spatial distribution of the data
#----------------------------------
plot.dat1 <- data.frame(coords_train)
plot.dat2 <- data.frame(coords_pred)
plot.dat1$status <- rep("Observed", nrow(plot.dat1))
plot.dat2$status <- rep("No Data", nrow(plot.dat2))
plot.dat <- rbind(plot.dat1, plot.dat2); names(plot.dat)
names(plot.dat) <- c("Longitude", "Latitude", "Status")
names(plot.dat)


data_locs <- ggplot(plot.dat, aes(x = Longitude, y = Latitude, color=Status)) +
  annotation_map(map_data("world"), colour = "light green", fill = "dark grey") +
  geom_point(size=1) +
  labs(color = "Status")+
  scale_color_manual(values=c('dark red', 'light green'))+
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        legend.position = "top", legend.spacing.x = unit(1.0, "cm"),
        legend.text = element_text(size = 16, color = "black"),
        legend.title = element_text(size = 16, face = "bold.italic"),
        text=element_text(size=16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)
  )
data_locs

##
rm(coords_pred, dat_pred) # free some memory spaces 
#-------Extract boundary for mesh
bnd <- inla.nonconvex.hull(as.matrix(coords),-0.03, -0.05, resolution=c(100,100))

###----Build non-convex hull mesh
mesh <- inla.mesh.2d(boundary = bnd, max.edge=c(0.6,4), cutoff=0.4)


# Visualise mesh
par(mfrow=c(1,1))
plot(mesh)  
points(coords, cex=0.1, col="red", pch=16)
mesh$n #----number of nodes

#----------------------------------------------------------------------------------
###---Build projector matrix A
A<-inla.spde.make.A(mesh=mesh,loc=as.matrix(coords));dim(A)


##---Create the SPDE
spde <- inla.spde2.matern(mesh, alpha=2)


##----specify the observation indices for estimation 
iset <- inla.spde.make.index(name = "spatial.field", spde$n.spde)

#-------Fitting the models----------------------------
#---All data
#--Recode for random effects 
datt$prov <- factor(datt$Prov_ID) # province random effects
datt$TYPE <- factor(datt$TYPE)
datt$set_typ <- as.factor(as.numeric(datt$TYPE))# settlement type random effects
# 1 - non-village; 2 - rural; 3 - urban
#table(datt$TYPE) 
#table(datt$set_typ) 


#---settlement - type and province nesting
Zsp <- as(model.matrix( ~ 0 + prov:set_typ, data = datt), "Matrix") 
datt$IDsp <- 1:nrow(datt)
datt$set_prov <- as.factor(apply(Zsp, 1, function(x){names(x)[x == 1]}))#--nesting


#--------------------------------------------------------------------------------
# prepare data for modelling
names(datt)
mod_covs <- 34:85 ; length(34:85)#----52 covariates

head(datt[, mod_covs]) #-----view the first sic rows of the covariates columns

# recode covariates (this is just for convenience)
colnames(datt)[mod_covs] <- paste0("x", 1:52, sep="") #--rename covariates 
names(datt)  #------check

# apply covariates scaling to the model covariates only
datt[,mod_covs] <- apply(datt[,mod_covs], 2, stdize)   

#---------------------------------------------------------------------------------
#      Fit building intensity model 
#---------------------------------------------------------------------------------
# select model covariates and key variables
best_covs_bld <- c("x3","x15","x27","x30", "x31","x33","x34","x35", "x36",
                   "x38","x46","x50","x52") # the best covariates selected based on stepwise algorithms


cov_bld <- datt[,c(best_covs_bld, "set_prov", "set_typ",  "prov", "IDsp")] # select only relevant variables

#---Build the stack for the bld
datt$BLDG21 <- datt$BLDG21 + 1# at least a value of 1 to allow for log transformation
stk_bld <- inla.stack(data=list(y=log(datt$BLDG21)), #the response
                      
                      A=list(A,1),  #the PROJECTION matrix
                      
                      effects=list(c(list(Intercept=1), #the Intercept
                                     iset),  #the spatial index
                                   #the covariates
                                   list(cov_bld)
                      ), 
                      #this is a quick name so you can call upon easily
                      tag='est_bld')


##--------------------------------------------------------------------------------------
# Model fitting (the best fit model)
# Model for building intensity (this is the best fit model following initial multiple iterations)
f_bld3 <- y ~ -1 + Intercept +  x3 + x15 + x27 + x30 + x31 + x33 + x34 + x35 + x36 + 
  x38 + x46 + x50 + x52 + f(spatial.field, model=spde) + f(IDsp, model='iid') + f(set_prov, model="iid")

#--------
# the function below takes about  4 minutes minutes to run on Intel(R) Core(TM) i5-6600 CPU @ 3.30 GHz; 32GB RAM Machine
system.time(bmod3<-inla(f_bld3, #the formula
                        data=inla.stack.data(stk_bld,spde=spde),  #the data stack
                        family= 'gaussian',   # probability distribution of the response variable
                        control.predictor=list(A=inla.stack.A(stk_bld),compute=TRUE),  #compute gives you the marginals of the linear predictor
                        control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                        verbose = FALSE)) #can include verbose=TRUE to see the log of the model runs
summary(bmod3) #----model summary

bind3 <-inla.stack.index(stk_bld, "est_bld")$data #--extract estimation indices 

# obtain the mean posterior estimates
bfit3 <- exp(bmod3$summary.linear.predictor[bind3,"mean"]) #--extract the backtransformed building intensity


#--extract fixed effects to 4 decimal places
(betab <- round(bmod3$summary.fix,4))


# Add the predicted building intensity to the data
datt$bld_pred <- bfit3

ls() # view contents of the workspace
rm(plot.dat, plot.dat1, plot.dat2, surv_df) # remove some irrelevant or duplicated files
#------------------------------------------------------------------------------------------------
#                      Fit population density model
#------------------------------------------------------------------------------------------------

#------------------------
# BHM  - uses the imperfectly observed settlement/building data directly 
#------------------------
datt$bld <- datt$BLDG21
datt$bld[datt$POPN==0] =NA # set all settlements without population to NA to allow for prediction
datt$dens_bhm <- datt$POPN/datt$bld
datt$dens_bhm[datt$dens_bhm==0] = 0.000001 # avoids mathematical issues with taking log of 0
hist(log(datt$dens_bhm))


# use the best covariates for 
best_covs_dens <- c("x3","x11","x15","x29","x30","x31","x32", "x35","x36","x39","x45","x46","x48","x51","x52")

####---Density Stack
covars_dens <- datt[,c(best_covs_dens,"set_prov", "set_typ","prov", "IDsp")]; dim(covars_dens)

#---Build the stack for the training set
##--(BHM)
stk_bhm<- inla.stack(data=list(y=datt$dens_bhm), #the response
                     
                     A=list(A,1),  #the PROJECTION matrix
                     
                     effects=list(c(list(Intercept=1), #the Intercept
                                    iset),  #the spatial index
                                  #the covariates
                                  list(covars_dens)
                     ), 
                     #this is a quick name so you can call upon easily
                     tag='est_bhm')

f_bhm <- y ~ -1 + Intercept +  x3 + x11 + x15 + x29 + x30 + x31 + x32 + x35 + x36 + 
  x39 + x45+ x46 + x48 + x51 + x52 + f(spatial.field, model=spde) + f(IDsp, model='iid') +
  f(set_typ, model='iid')

# the function below takes about 8 minutes to run on Intel(R) Core(TM) i5-6600 CPU @ 3.30 GHz; 32GB RAM Machine
system.time(mod_bhm <-inla(f_bhm, #the formula
                           data=inla.stack.data(stk_bhm,spde=spde),  #the data stack
                           family= 'gamma',   #which family the data comes from
                           control.predictor=list(A=inla.stack.A(stk_bhm),compute=TRUE),  #compute gives you the marginals of the linear predictor
                           control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                           verbose = FALSE)) #can include verbose=TRUE to see the log of the model runs
#summary(mod_bhm) #----model summary

ind_bhm <-inla.stack.index(stk_bhm, "est_bhm")$data #--estimation indices 
fit_bhm <- exp(mod_bhm$summary.linear.predictor[ind_bhm,"mean"]) #--extract the backtransformed pop_hat

sum(fitbhm <- fit_bhm*datt$bld, na.rm=T) # check total count

#--------------------------------------------------
# TSBHM - uses the bias-adjusted settlement/building data
#-------------------------------------------------------
datt$dens_tsbhm <- datt$POPN/datt$bld_pred # uses the predicted (bias-corrected) building intensity
datt$dens_tsbhm[datt$dens_tsbhm==0] =  0.000001 # avoids mathematical issues with taking log of 0

# use the best covariates for 
best_covs_dens <- c("x3","x11","x15","x29","x30","x31","x32", "x35","x36","x39","x45","x46","x48","x51","x52")
####---Density Stack
covars_dens <- datt[,c(best_covs_dens,"set_prov", "set_typ","prov", "IDsp")]; dim(covars_dens)

#---Build the stack for the training set
##--(tsbhm)
stk_tsbhm<- inla.stack(data=list(y=datt$dens_tsbhm), #the response
                       
                       A=list(A,1),  #the PROJECTION matrix
                       
                       effects=list(c(list(Intercept=1), #the Intercept
                                      iset),  #the spatial index
                                    #the covariates
                                    list(covars_dens)
                       ), 
                       #this is a quick name so you can call upon easily
                       tag='est_tsbhm')


f_tsbhm <- y ~ -1 + Intercept +  x3 + x11 + x15 + x29 + x30 + x31 + x32 + x35 + x36 + 
  x39 + x45+ x46 + x48 + x51 + x52 + f(spatial.field, model=spde) + f(IDsp, model='iid') +
  f(set_typ, model='iid')

# the function below takes approximately 7 minutes to run on Intel(R) Core(TM) i5-6600 CPU @ 3.30 GHz; 32GB RAM Machine
system.time(mod_tsbhm <-inla(f_tsbhm, #the formula
                             data=inla.stack.data(stk_tsbhm,spde=spde),  #the data stack
                             family= 'gamma',   #which family the data comes from
                             control.predictor=list(A=inla.stack.A(stk_tsbhm),compute=TRUE),  #compute gives you the marginals of the linear predictor
                             control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                             verbose = FALSE)) #can include verbose=TRUE to see the log of the model runs

summary(mod_tsbhm) #----model summary
 
ind_tsbhm <-inla.stack.index(stk_tsbhm, "est_tsbhm")$data #--estimation indices 
fit_tsbhm <- exp(mod_tsbhm$summary.linear.predictor[ind_tsbhm,"mean"]) #--extract the backtransformed pop_hat
sum(fittsbhm <- fit_tsbhm*datt$bld_pred, na.rm=T) 


##--------------------------------------------------------------------------------------------
#####-----POSTERIOR SIMULATION-----------------------------------------------------
#------------------------------------------------------------------------------------------
# extract posterior parameter estimates of random effects
stypp <- function(mod, dat)
{
  dat$stype <- rep(1, nrow(dat))
  st <- mod$summary.random$set_typ$mean
  uniq <- unique(dat$set_typ)
  uniq[1]
  for(i in  1:nrow(dat))
  {
    
    for(j in 1:3)
    {
      if(dat$set_typ[i]==uniq[j]) dat$stype[i] = st[j]
    }
  }
  dat$stype
}

# duplicate data
datt_bhm <- datt

# settlement type random effects
datt_bhm$set_re <- stypp(mod_bhm, datt_bhm)
datt$set_re <- stypp(mod_tsbhm, datt)

# bhm
datt_bhm$BLD <- datt_bhm$bld

# tsbhm
datt$BLD <- datt$bld_pred



# function for posterior simulation and prediction-------------
#-----------------------------------------------------------
# Draws parameter samples from the posterior in each iteration
# Uses each set of drawn parameter values to make predictions of the response
# Obatins and stores aggregated totals in each iteration 
#----------------------------------------------------------

simPops_gg <- function(model, dat, Aprediction, run)
{
  
  #--------------------------------------------------
  # model: the inla model from where samples will be drawn 
  # dat: the data frame containing the prediction covariates
  # Aprediction: projection matrix building over the prediction domains
  # run: number of iterations (any value from 30 and above may suffice)
  # No need to worry about convergence as INLA is deterministic and samples from the stationary distribution
  #----------------------------------------
 
  fixedeff  <- dens_hat <- pop_hat <- matrix(0, nrow=nrow(dat), ncol = run)
  #inla.seed = as.integer(runif(1)*.Machine$integer.max) ## generates sample seed
  inla.seed =  481561959 # simulation seed for reproducible samples
  set.seed(inla.seed)
  print(inla.seed)
  m1.samp <- inla.posterior.sample(run, model, seed = inla.seed ,selection=list(x3=1, x11=1, x15=1, # initialising the coefficients
                                                                                x29=1, x30=1, x31=1,
                                                                                x32=1, x35=1, x36=1,
                                                                                x39=1, x45=1, x46=1,
                                                                                x48=1, x51=1, x52=1),num.threads="1:1")
  
  
  sfield_nodes_mean <- model$summary.random$spatial.field['mean'] # extract the spatial random effect of the model
  field_mean <- (Aprediction%*% as.data.frame(sfield_nodes_mean)[, 1])
  for(i in 1:run)
  {
    
    fixedeff[,i] <-  
      model$summary.fixed['Intercept', 'mean'] + # allow a fixed intercept (this can also vary)
      m1.samp[[i]]$latent[1,] * dat[,'x3'] +
      m1.samp[[i]]$latent[2,] * dat[,'x11'] +
      m1.samp[[i]]$latent[3,] * dat[,'x15'] +
      m1.samp[[i]]$latent[4,] * dat[,'x29'] +
      m1.samp[[i]]$latent[5,] * dat[,'x30'] + 
      m1.samp[[i]]$latent[6,] * dat[,'x31'] +
      m1.samp[[i]]$latent[7,] * dat[,'x32'] +
      m1.samp[[i]]$latent[8,] * dat[,'x35'] +
      m1.samp[[i]]$latent[9,] * dat[,'x36'] +
      m1.samp[[i]]$latent[10,] * dat[,'x39'] +
      m1.samp[[i]]$latent[11,] * dat[,'x45'] +
      m1.samp[[i]]$latent[12,] * dat[,'x46'] +
      m1.samp[[i]]$latent[13,] * dat[,'x48'] +
      m1.samp[[i]]$latent[14,] * dat[,'x51'] +
      m1.samp[[i]]$latent[15,] * dat[,'x52'] +
      
      
      model$summary.random$IDsp['mean'][,1] + # iid random effect
      dat$set_re + # settlement type random effect
      
      field_mean[,1] # spatial random effect
    
    dens_hat[,i]<- exp(fixedeff[,i]) # predicted population density
    pop_hat[,i] <- dens_hat[,i]*dat$BLD # predicted population count
  }
  
  # Obtain posterior inference
  mean_dens_hat <- apply(dens_hat, 1, mean, na.rm=T) # mean density
  mean_pop_hat <- apply(pop_hat, 1, mean, na.rm=T) # mean population count 
  median_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.5), na.rm=T) # median population count 
  lower_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.025), na.rm=T) # lower bound of population count (95% credible interval)
  upper_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.975), na.rm=T) # upper bound of population count (95% credible interval)
  sd_pop_hat <- apply(pop_hat, 1, sd, na.rm=T) # standard deviation of population count (95% credible interval)
  uncert_pop_hat <- (upper_pop_hat - lower_pop_hat)/mean_pop_hat# explicit estimates of uncertainties of estimated population counts 
  
  
  # add to the dataset
  dat$mean_dens_hat <- mean_dens_hat
  dat$mean_pop_hat <- mean_pop_hat
  dat$median_pop_hat <- median_pop_hat
  dat$lower_pop_hat <- lower_pop_hat
  dat$upper_pop_hat <- upper_pop_hat
  dat$uncert_pop_hat <- uncert_pop_hat
  dat$sd_pop_hat <- sd_pop_hat
  
  
  # save output as a list
  output <- list(pop_hat = pop_hat, # this is an nrow(dat) by run dimension matrix of posterior samples
                 est_data = dat)
  
}

run=2000 # number of iterations 
system.time(str(sim.bhm <- simPops_gg(mod_bhm,datt_bhm, A, run)))  # BHM - takes approximately 6 minutes 

system.time(str(sim.tsbhm <- simPops_gg(mod_tsbhm,datt, A, run))) # TSBHM- takes approximately 6 minutes 

#2065842720;1128763477
# Visualise the traceplots of the posterior samples 
# Load ggplot2
library(ggplot2)

##-----Posterior samples traceplots of 6 randomly selected census units
samp <- sample(nrow(datt), 6)
par(mfrow=c(3,2), mar=c(5,5,2,1))
for(j in samp)
{
  plot.ts(sim.tsbhm$pop_hat[j,], ylab = "pop_hat", col="steel blue",
          cex.main=2,  main=paste0("Pop_hat samples for CU ", datt$CU_Name[j], sep=""))
  abline(h=mean(sim.tsbhm$pop_hat[j,]), lwd=2, col=2)
}



#------------------------------------------------------------
###---CU level estimates
#----------------------------------------------------------

# check national totals
sum(sim.bhm$est_data$mean_pop_hat, na.rm=T)    # bhm
sum(sim.tsbhm$est_data$mean_pop_hat, na.rm=T)  # tsbhm

# Exract CU level posterior estimates as data frames with other variables
cu.tsbhm <- sim.tsbhm$est_data # TSBHM 
cu.bhm <- sim.bhm$est_data # BHM 

#-----------------------------------
#   Calculate  Relative Error rate 
#-------------------------------
# ----TSBHM
predicted.total.tsbhm <- sum(cu.tsbhm$mean_pop_hat[!is.na(cu.tsbhm$POPN)], na.rm=T) # predicted total
observed.total <- sum(cu.tsbhm$POPN,na.rm=T) # observed total 

# error rate
erate.tsbhm <- abs(predicted.total.tsbhm-observed.total)/observed.total


# ---- BHM
predicted.total.bhm <- sum(cu.bhm$mean_pop_hat[!is.na(cu.bhm$POPN)], na.rm=T) # predicted 
observed.total <- sum(cu.bhm$POPN,na.rm=T)

# error rate
erate.bhm <- abs(predicted.total.bhm-observed.total)/observed.total


# --- relative error rate
rerate <- erate.tsbhm/erate.bhm

## --- reduction in relative error rate (%)
(1-rerate)*100 # ~32%
#-------------------------------------------------------------------------------------------------


#####----Join the posterior draws to the dataframe------------------------------------------------------
dim(datt)
names(datt)
data.all <- datt[,c("Prov_Name" ,"Dist_Name", "LLG_Name")]

# TSBHM
data.tsbhm <- data.frame(cbind(data.all, sim.tsbhm$pop_hat))#---Contains all the simulated posterior matrix
dim(data.tsbhm);names(data.tsbhm)


# BHM
data.bhm<- data.frame(cbind(data.all, sim.bhm$pop_hat))#---Contains all the simulated posterior matrix
dim(data.bhm);names(data.bhm)


#-----------------------------------------------------------------------------
#   Obtain Admin totals with uncertainties
#----------------------------------------------------------------------------

thin=5 #----Thinning choses every 5th sample for posterior inference 
thinn <- seq((4+0.2*run),run, by=thin) #--There is a burnin period of the first 20% of the total sample

#----NATIONAL------------------------------
nat_total <- function(dat, thinn)
{
  p_hat <- dat[,thinn]
  tots <- apply(p_hat,2, sum, na.rm=T) #Col sums
  
  tot_sd  <- sd(tots, na.rm=T)
  
  tot_mean  <- mean(tots, na.rm=T)
  
  tot_lower <- quantile(tots, probs=c(0.025))
  tot_median <- quantile(tots, probs=c(0.5))
  tot_upper <- quantile(tots, probs=c(0.975))
  
  return(estimates <- data.frame(estimates=unlist(list(total=tot_mean, 
                                                       lower=tot_lower, 
                                                       median=tot_median, 
                                                       upper=tot_upper))))
}
(national.tsbhm <- nat_total(data.tsbhm, thinn)) # TSBHM
#(national.bhm<- nat_total(data.bhm, thinn)) # BHM



#write.csv(national1, file=paste0(results_path, "/updated2/National_estimates_main_gamma-gaussian.csv"))

#  Province total and uncertainties
prov_est_gg <- function(datp, thinn)
{
  provnames <- unique(datp$Prov_Name)
  outP <- matrix(0, nrow=length(provnames), ncol=3)
  for(j in 1:length(provnames))
  {
    prov <- datp[datp$Prov_Name==provnames[j],]
    #ptots <- apply(prov[,6:(5+run)], 2, sum, na.rm=T)
    ptots <- apply(prov[,thinn], 2, sum, na.rm=T)
    ptots_sd <- sd(ptots, na.rm=T)
    
    ptot_mean1  <- mean(ptots, na.rm=T)
    ptot_lower <- quantile(ptots, probs=c(0.025))
    #ptot_meadian <- quantile(ptots, probs=c(0.5))
    ptot_upper <- quantile(ptots, probs=c(0.975))
    
    pestimates <- round(c(ptot_mean1, ptot_lower, ptot_upper), 2)
    outP[j,] <- pestimates
  }
  outP <- data.frame(outP)
  return(PROV_est <- data.frame(names = provnames,
                                total = outP[,1],
                                lower = outP[,2],
                                #median = outP[,3],
                                upper = outP[,3]))
}

prov.tsbhm <- prov_est_gg(data.tsbhm,thinn) # TSBHM
prov.bhm <- prov_est_gg(data.bhm,thinn) 

#------------------------------------------------------------
###---District level estimates 
#----------------------------------------------------------
dist_est_gg <- function(datd, thinn)
{
  dnames <- unique(datd$Dist_Name)
  outd <- matrix(0, nrow=length(dnames), ncol=3)
  
  for(j in 1:length(dnames))
  {
    dist <- datd[datd$Dist_Name==dnames[j],]
    dtots <- apply(dist[,thinn], 2, sum, na.rm=T)
    
    dtots_sd <- sd(dtots, na.rm=T)
    
    dtot_mean1  <- mean(dtots, na.rm=T)
    
    dtot_lower <- quantile(dtots, probs=c(0.025))
    dtot_upper <- quantile(dtots, probs=c(0.975))
    destimates <- round(c(dtot_mean1, dtot_lower, dtot_upper), 2)
    outd[j,] <- destimates
  }
  outd <- data.frame(outd)
  return(dist_est <- data.frame(names = dnames,
                                total = outd[,1],
                                lower = outd[,2],
                                upper = outd[,3]))
}

(dist1 <- dist_est_gg(data.tsbhm,thinn)) # TSBHM

#write.csv(dist1, file=paste0(res_path, "/District_data_TSBHM.csv"), row.names=F)
#------------------------------------------------------------
###---LLG level estimates 
#----------------------------------------------------------
llg_est_gg <- function(datl, thinn)
{
  lnames <- unique(datl$LLG_Name)
  outl <- matrix(0, nrow=length(lnames), ncol=3)
  
  for(j in 1:length(lnames))
  {
    llg <- datl[datl$LLG_Name==lnames[j],]
    ltots <- apply(llg[,thinn], 2, sum, na.rm=T)
    ltots_sd <- sd(ltots, na.rm=T)
    
    ltot_mean1  <- mean(ltots, na.rm=T)
    ltot_lower <- quantile(ltots, probs=c(0.025))
    ltot_upper <- quantile(ltots, probs=c(0.975))
    
    lestimates <- round(c(ltot_mean1, ltot_lower, ltot_upper), 3)
    outl[j,] <- lestimates
  }
  outl <- data.frame(outl)
  return(llg_est <- data.frame(names = lnames,
                               total = outl[,1],
                               lower = outl[,2],
                               upper = outl[,3]))
}

llg1 <- llg_est_gg(data.tsbhm,  thinn) # TSBHM 


library(ggpubr)

##
var2Add <- c("Prov_Name", "Dist_Name", "LLG_Name", "Ward_Name", "CU_Name", "POPN", "BLDG21", "bld_pred",
             "set_typ", "mean_dens_hat", "mean_pop_hat", "lower_pop_hat", "upper_pop_hat", "uncert_pop_hat")

# scatter and 2d plots
## TSBHM 


# TSBHM
dim(dat_tsbhm <- cu.tsbhm) 
dat_tsbhm$method <- rep("TSBHM", nrow(dat_tsbhm))
dat_tsbhm$ldens <- log(dat_tsbhm$mean_dens_hat) # log posterior pop density


# BHM
dim(dat_bhm <-cu.bhm) 
dat_bhm$method <- rep("BHM", nrow(dat_bhm))
dat_bhm$ldens <- log(dat_bhm$mean_dens_hat)  # log posterior pop density


#-----------------------------------------------------------------------------
###  FIGURE 5 (MAIN MANUSCRIPT)
#----------------------------------------------------------------------------------
plot_cu1 <- ggplot(cu.tsbhm, aes(POPN, mean_pop_hat)) +
  geom_pointrange(aes(ymin = lower_pop_hat, ymax = upper_pop_hat), size=0.5)+
  # geom_errorbar(aes(ymin = lower, ymax = upper), size=1,width = 0.5, colour='dark grey')+
  geom_line(size=1, colour='red') +
  theme_bw()+
  #facet_wrap(~method, scales="free", nrow=1)+
  theme(strip.text = element_text(size = 15),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14))#+

plot_cu1<- ggpar(plot_cu1, ylab="Predicted Population", xlab="Observed Population",
                 legend = "right", legend.title=element_text("Density",size=22) ,
                 font.legend=c(20),
                 font.label = list(size = 20, face = "bold", color ="red"),
                 font.x = c(20),
                 font.y = c(20),
                 font.main=c(16),
                 font.xtickslab =c(20),
                 font.ytickslab =c(20),
                 #xticks.by = TRUE,
                 #xlim= c(0, max(df2$x)),
                 xtickslab.rt = 45, ytickslab.rt = 45)

plot_cu1
#''''''


##  district summary-----------------
(dst.data1 <- cu.tsbhm %>% group_by(Dist_Name) %>%
   summarise(mean = sum(mean_pop_hat, na.rm=T),
             pop = sum(POPN, na.rm=T)))
dst.data1$method <- rep("TSBHM", nrow(dst.data1))

##########

dstt <- cu.tsbhm %>% group_by(Dist_Name)%>% drop_na(mean_pop_hat) %>%
  summarise(Obs = sum(POPN, na.rm=T),
            Pred = round(sum(mean_pop_hat, na.rm=T)))

cup1 <-ggplot(data=dstt,aes(Obs,Pred)) + # swapped 
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black',
                 contour_var = "ndensity") + 
  scale_fill_continuous(type="viridis") +
  geom_smooth(method=lm,linetype=2,colour="red",se=F) + 
  scale_y_continuous(breaks=seq(0, max(dstt$Obs),50000))+
  guides(alpha="none") +
  geom_point() + commonTheme

plot_dst1<- ggpar(cup1, xlab="Observed Population", ylab="Predicted Population",
                  legend = "right", legend.title=element_text("Density",size=22) ,
                  font.legend=c(20),
                  legend.show=F,
                  font.label = list(size = 20, face = "bold", color ="red"),
                  font.x = c(20),
                  font.y = c(20),
                  font.main=c(16),
                  font.xtickslab =c(20),
                  font.ytickslab =c(20),
                  #xticks.by = TRUE,
                  xlim= c(0, max(dstt$Obs)),
                  xtickslab.rt = 45, ytickslab.rt = 45)

plot_dst1


ggarrange(plot_cu1, plot_dst1+ rremove("x.text"), 
          labels = c("(A)", "(B)"),
          nrow=2)


#----------------------------------------------------------------------------------------------
#             FIGRE U -----------------------------------------------------------------------
# barplots
prov1$method <- rep("TSBHS", nrow(prov1))
prov2$method <- rep("BHS", nrow(prov2))



prv_dat <- bind_rows(prov1, prov2)
pcp <- ggplot(prv_dat, aes(x=names, y=total, group=method, fill=method)) +
  #geom_line(linetype="dashed") +
  geom_bar(aes(group=method), alpha=0.5, stat="identity") +
  geom_pointrange( aes(ymin=lower, ymax=upper, color=method), alpha=0.9, size=0.8)+
  scale_fill_manual(values=c("steel blue", "light blue"))

pdcp <- ggpar(pcp, ylab="Predicted Total Population", xlab="Province",
              legend = "right", legend.title = "Method",
              font.label = list(size = 18, face = "bold", color ="red"),
              font.x = c(16),
              font.y = c(16),
              xtickslab.rt = 45, ytickslab.rt = 45)

pdcp

##------------------- MORE CSATTER PLOTS -----------------------------------------------
dst <- ggscatter(dst.data1, x = "pop", y = "mean",
                 color = "black", shape = 16, size = 2, # Points color, shape and size
                 add = "reg.line",  # Add regressin line
                 #ellipse = TRUE,
                 add.params = list(color = "red", fill = "light green"), # Customize reg. line
                 conf.int = TRUE, # Add confidence interval
                 cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                 cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)
dst11 <- ggpar(dst, xlab="Observed population count", ylab="Predicted population count",
               legend = "right", legend.title=element_text("Method",size=22),
               font.legend=c(18),
               palette = c("Blue"),
               font.label = list(size = 18, face = "bold", color ="red"),
               font.x = c(22),
               font.y = c(18),
               font.main=c(22),
               font.xtickslab =c(18),
               font.ytickslab =c(16),
               #ylim=c(0, 13000),
               xtickslab.rt = 45, ytickslab.rt = 45)
dst11




prv <- ggscatter(prv.data1, x = "pop", y = "mean",
                 color = "black", shape = 16, size = 2, # Points color, shape and size
                 add = "reg.line",  # Add regressin line
                 #ellipse = TRUE,
                 add.params = list(color = "red", fill = "steel blue"), # Customize reg. line
                 conf.int = TRUE, # Add confidence interval
                 cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                 cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)

prv11 <- ggpar(prv, xlab="Observed population count", ylab="Predicted population count",
               legend = "right", legend.title=element_text("Method",size=22),
               font.legend=c(18),
               palette = c("red"),
               font.label = list(size = 18, face = "bold", color ="red"),
               font.x = c(22),
               font.y = c(18),
               font.main=c(22),
               font.xtickslab =c(18),
               font.ytickslab =c(16),
               #xticks.by = 165319,
               #xlim=c(min(prv.data1$pop), max(prv.data1$pop)),
               xtickslab.rt = 45, ytickslab.rt = 45)
prv11



#################################################################################
#####-------------------Cross validation
####################################################################################

#----Extract settement type effects
set_t <- function(dat, st)
{
  #st <- mod$summary.random$set_typ$mean
  uniq <- unique(dat$set_typ)
  uniq[1]
  for(i in  1:nrow(dat))
  {
    
    for(j in 1:3)
    {
      if(dat$set_typ[i]==uniq[j]) dat$set_typ2[i] = st[j]
    }
    
  }
  dat$set_typ2
}

#------------------------------------------------

# In-Sample 

#-------------------------------------------------

# #---------------------------------------------
#   IN-SAMPLE CROSS VALIDATION
# #-----------------------------------------

# Summarise posterior estimates for cross-validation

## bhm
bhm_post <- data.frame(obs = cu.bhm$POPN, # observed count
                       pred = cu.bhm$mean_pop_hat, # predicted count
                       dens = cu.bhm$dens_bhm, # observed density
                       densp = cu.bhm$mean_dens_hat) %>%# predicted population density
drop_na() # remove all NAs

## tsbhm
tsbhm_post <- data.frame(obs = cu.tsbhm$POPN, # observed count
                       pred = cu.tsbhm$mean_pop_hat, # predicted count
                       dens = cu.tsbhm$dens_tsbhm, # observed density
                       densp = cu.tsbhm$mean_dens_hat) %>%# predicted population density
  drop_na() # remove all NAs


### count
(met11 <- mod_metrics2(bhm_post$obs, 
                       bhm_post$pred)) #bhm

(met22 <- mod_metrics2(tsbhm_post$obs, 
                       tsbhm_post$pred)) #tsbhm

(met_pop <- t(data.frame(bhm = unlist(met11),
                         tsbhm = unlist(met22)))) # join

### densities 
(met11d <- mod_metrics2(bhm_post$dens, 
                        bhm_post$densp)) # bhm

(met12d <- mod_metrics2(tsbhm_post$dens, 
                        tsbhm_post$densp))# tsbhm


(met_dens <- t(data.frame(bhm = unlist(met11d),
                         tsbhm = unlist(met12d)))) #join

####----------------------------------------------------------------------------
##
#   OUT-OF-SAMPLE K-FOLD CROSS-VALIDATION 
#
#--------------------------------------------------------
set.seed(9444506) # set sampling seed                                    

dat = datt
k_folds = 5
N <- nrow(dat)
######
ind_train <- factor(sample(x = rep(1:k_folds, each = floor(N/ k_folds)),  # Sample IDs for training data
                           size = N, rep=T))

table(as.numeric(ind_train)) 
dat$k_fold <- as.numeric(ind_train)
coords <- cbind(dat$lon, dat$lat)


k_uniq <-sort(unique(dat$k_fold))
#metrics_cv <- matrix(0, nrow=length(k_uniq), ncol=5)#5 metrics for each fold

met_list_out <- list()
for(i in 1:length(k_uniq))
{
  
  print(paste0("fold_", i, sep=""))
  train_ind <- which(dat$k_fold!=k_uniq[i])
  dim(train <- dat[train_ind, ])#---train set for fold i
  dim(test <- dat[-train_ind, ]) #---test set for fold i
  
  train_coords <- coords[train_ind,]
  test_coords <- coords[-train_ind,]
  
  
  ###---Create projection matrices for training and testing datasets
  Ae<-inla.spde.make.A(mesh=mesh,loc=as.matrix(train_coords));dim(Ae) #training
  
  
  
  ########################
  covars_train <- train[,c("x3","x11","x15","x29","x30", "x31","x32","x35", "x36",
                           "x39","x45","x46","x48","x51","x52","set_prov", "set_typ", 
                           "prov", "IDsp")]; dim(covars_train)
  
  train$dens2 <- train$dens_tsbhm
  #---Build the stack for the training set
  stk_train <- inla.stack(data=list(y=train$dens2), #the response
                          
                          A=list(Ae,1),  #the A matrix; the 1 is included to make the list(covariates)
                          
                          effects=list(c(list(Intercept=1), #the Intercept
                                         iset),  #the spatial index
                                       #the covariates
                                       list(covars_train)
                          ), 
                          #this is a quick name so you can call upon easily
                          tag='train')
  
  
  
  covars_test <- test[,c("x3","x11","x15","x29","x30", "x31","x32","x35", "x36",
                         "x39","x45","x46","x48","x51","x52","set_prov", "set_typ", 
                         "prov", "IDsp")]; dim(covars_test)

  
  ###---Rerun INLA for model test prediction
  model <-inla(f_tsbhm, #the formula
               data=inla.stack.data(stk_train,spde=spde),  #the data stack
               family= 'gamma',   #which family the data comes from
               control.predictor=list(A=inla.stack.A(stk_train),compute=TRUE),  #compute gives you the marginals of the linear predictor
               control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
               verbose = FALSE) #can include verbose=TRUE to see the log of the model runs
  summary(model)
  
  
  
  mod = mod_tsbhm
  #######################----Extract Spatial Random effects
  sfield_nodes_mean <- mod$summary.random$spatial.field['mean']
  field_mean <- (A%*% as.data.frame(sfield_nodes_mean)[, 1])
  
  
  
  
  ###----Extract settlement type random effects
  set <- set_t(dat, mod$summary.random$set_typ$mean)
  
  ##--------
  fixed <-  
    model$summary.fixed['Intercept', 'mean'] +
    model$summary.fixed['x3', 'mean']* test[,'x3'] +
    model$summary.fixed['x11', 'mean'] * test[,'x11'] +
    model$summary.fixed['x15', 'mean'] * test[,'x15'] +
    model$summary.fixed['x29', 'mean'] * test[,'x29'] +
    model$summary.fixed['x30', 'mean'] * test[,'x30'] + 
    model$summary.fixed['x31', 'mean'] * test[,'x31'] +
    model$summary.fixed['x32', 'mean'] * test[,'x32'] +
    model$summary.fixed['x35', 'mean'] * test[,'x35'] +
    model$summary.fixed['x36', 'mean'] * test[,'x36'] +
    model$summary.fixed['x39', 'mean'] * test[,'x39'] +
    model$summary.fixed['x45', 'mean'] * test[,'x45'] +
    model$summary.fixed['x46', 'mean'] * test[,'x46'] +
    model$summary.fixed['x48', 'mean',] * test[,'x48'] +
    model$summary.fixed['x51', 'mean'] * test[,'x51'] +
    model$summary.fixed['x52', 'mean'] * test[,'x52'] +
    
    mod$summary.random$IDsp['mean'][-train_ind,1] +
    set[-train_ind] +
    field_mean[-train_ind,1]
  
  dens_ht <- exp(fixed)
  sum(pop_ht <- dens_ht*test$bld_pred)
  
  
  ### scatter plots 
  par(mfrow=c(1,1))
  plot(test$dens_tsbhm, dens_ht, xlab = "Observed", 
       ylab = "Predicted", col=c('dark green','orange'),
       pch=c(16,16), cex.axis=1.5)
  abline(0,1)
  legend("topleft", c("Observed", "Predicted"), col=c("dark green", "orange"), pch=c(16,16),
         bty="n", cex=1.5) 
  
  # calculate fit metrics
  (met_out <- mod_metrics2(test$dens_tsbhm,  
                           dens_ht))
  
  met_list_out[[i]]<- unlist(met_out)
}
met_list_out_dat <- do.call(rbind,met_list_out)
metrics_out <- apply(met_list_out_dat, 2, mean)


# compare results 
#rbind(met_dens[2,], metrics_out)
##----Save workspace
#save.image(paste0(results_path, "/workspace_gg_model.Rdata"))
