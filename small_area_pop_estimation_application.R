####--TITLE:Small area population estimation from health intervention campaign surveys and partially observed settlement data--#####
####--METHODS: GEOSTATISTICAL BAYESIAN HIERARCHICAL REGRESSION MODEL
####--AUTHOR: DR CHIBUZOR CHRISTOPHER NNANATU
####--INSTITUTION: WORLDPOP, UNIVERSITY OF SOUTHAMPTON 
####--DATE:19 DECEMBER 2022. REVIEWED: 25/03/2024
###=====================================================================================================

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
###--Specify various file paths (set input and output paths)
#--------------------------------------------------------------------------------
path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP000008_UNFPA_PNG/Working/Chris/paper/paper1/real_data/"
surv_path <- paste0(path, "input_data/survey/")
shape_path <- paste0(path, "input_data/boundary/")
results_path <- paste0(path, "results/")

#-------------------------------------------------------------------------------
#------------------------Load datasets-----------------------------------------        
#--------------------------------------------------------------------------------
shp_png <- st_read(paste0(shape_path,"CU/PNG_CU_32100_B.shp")) #load the CU shapefile

load(paste0(surv_path,"ModelDataPNG_BI.RData")) # load the .Rdata containing the survey data


#-----------------------------------------------------------------------------
#--------------------------------------------------------------------------
###---Explore  and clean data ----------------------------------------
#--------------------------------------------------------------------------
ls() # check the contents of the file
    # the survey dataframe is called 'covs' 


names(covs); str(covs)
covDF <- as.data.frame(covs) # data should be a data frame
#plot(covDF$lon, covDF$lat)
#names(covDF); str(covDF); dim(covDF); dim(shp)# shp has one extra row

# clean 
shp_png <- shp_png[-32100,] #--remove row 32100 from the shapefile as it does not exist
#plot(shp["Shape_Area"])
dim(covDF); dim(shp_png) #---confirm ir is now same number of rows as the survey data

# convert sf shapefile to sp 
shp <- as(st_geometry(shp_png), "Spatial") # converts sf data to sp. 
plot(shp)


#--------------------------------------------------------------------------
##----Extract the coordiates - the centroids of the CUs
#--------------------------------------------------------------------
library(sp)
covDF$lon <- coordinates(shp)[,1]#--add lon-lat to the data
covDF$lat <- coordinates(shp)[,2]
covDF2 <- covDF

#---------------------------------------------------------------------------
# covariates Z-score scaling 
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
  WAIC = mod$waic$waic
  DIC = mod$dic$dic
  CPO = -sum(log(mod$cpo$cpo), na.rm=T)
  output <- list(MAE  = INACCURACY ,
                 RMSE = RMSE,
                 BIAS = abs(BIAS),
                 COR=COR)
  return(output)
}

#--------------------------------------------------------------------------
# More data preparation s
#----------------------------------------------------------------------------
datt <- covDF2 #-------rename data
mod_covs <- 34:85 ; length(34:85)#----COVARIATES COLUMNS in datt

head(datt[, mod_covs]) #-----view the covariates columns

# apply covariates scaling to only the model covariates
datt[,mod_covs] <- apply(datt[,mod_covs], 2, stdize)   

# explore the datasets 
par(mfrow=c(1,1))
hist(datt$POPN)  #---observed pop total per CU
boxplot(datt$POPN) 
hist(datt$BLDG21) #---Building intensity ; hist(log(datt$BLDG21))
#------------------------------------------------------------------------

# run covariates selection 
colnames(datt)[mod_covs] <- paste0("x", 1:52, sep="") #--rename covariates 
names(datt)  #------confirm

# divide data into observed/training and unobserved/test
dim(dat_train <- datt[!is.na(datt$POPN),]) #-- Training set - all CUs with Observations only (16872 CUs)
dim(dat_pred <- datt[is.na(datt$POPN),])   #-- Prediction set - all unsampled CUs (15095 CUs)

#----------------Covariates Selection-----------------------------------
# Run GLM-based stepwise selection 

#install.packages("car")
library(car) ##--For calculating variance inflation factor (vif)
library(dplyr)
library(tidyverse)


#######-------For Density================================================================
names(dat_train)
covs_pop <- dat_train[,c(91,mod_covs)] #--subset for variables selection
covs_pop1 <- covs_pop %>% drop_na() #- model covariates only without NAs
dim(covs_pop1 <- covs_pop1[is.finite(covs_pop1$POPN),]) #--checks

#-----------------FITTIG THE STEPWISE REG MODELS-------------------------------------
fpop <- glm.nb(POPN ~., data=covs_pop1) #---negative binomial

step_pop <- stepAIC(fpop, scale = 0,
                    direction = c("both"),
                    trace = 1, keep = NULL, steps = 1000, use.start = FALSE,
                    k = 2)
step_pop

# only sttaistocally significant variables with vif less than 5 are retained
vif_pop = vif(fpop2)
vif_pop[which(vif_pop < 5)]

# Best fit covariates eventually selected
cov_pop <- c("x3","x11","x15","x29","x30", "x31","x32","x35", "x36",
             "x39","x45","x46","x48","x51","x52")

cov_names <- names(covDF2)[mod_covs]
covariates_pop <- cov_names[c(3, 11, 15, 29, 30, 31, 32, 
                              35, 36, 39, 45, 46, 48, 51, 52)]
pop_covs <- data.frame(cov = cov_pop, name=covariates_pop)
#-------------------------------------

##----Make correlation plot
#install.packages("corrplot")
require(corrplot)
png(paste0(results_path,"cor_plots.png"))
corrplot(
  cor(covs_pop[,cov_pop]),
  method = 'square',
  type = 'upper',
  tl.col = 'black',
  tl.cex = 1,
  col = colorRampPalette(c('purple', 'dark green'))(200)
)
dev.off()


#---------------------------------------------------------------------
#           BAYESIAN STATISTICAL MODELLING
#-------------------------------------------------------------------

# Define rhe points for both training and testint/predictions separately
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


###----Build non-hull mesh
bnd <- inla.nonconvex.hull(as.matrix(coords),-0.03, -0.05, resolution=c(100,100))
mesh <- inla.mesh.2d(boundary = bnd, max.edge=c(0.6,4), cutoff=0.4)


# Visualise mesh
par(mfrow=c(1,1))
#png(paste0(results_path, "/plots/mesh.png"))
plot(mesh)  
points(coords, cex=0.1, col="red", pch=16)
#dev.off()
mesh$n #----number of nodes

#----------------------------------------------------------------------------------
###---Build projector matrix A
A<-inla.spde.make.A(mesh=mesh,loc=as.matrix(coords));dim(A)


##---Create the SPDE
spde <- inla.spde2.matern(mesh, alpha=2)


##----specify the observation indices for estimation 
iset <- inla.spde.make.index(name = "spatial.field", spde$n.spde)

#-------Fitting the models----------------------------

#--Recode for random effects 
#---All data
datt$prov <- datt$Prov_ID # province random effects
datt$set_typ <- as.factor(as.numeric(datt$TYPE)) # settlement type random effects


# Add settlement type and province interactions
Zsp <- as(model.matrix( ~ 0 + prov:set_typ, data = datt), "Matrix") 
datt$IDsp <- 1:nrow(datt)
datt$set_prov <- as.factor(apply(Zsp, 1, function(x){names(x)[x == 1]}))#--nesting

#---------------------------------------------------------------------------------
#      Gaussian-Gamma
#---------------------------------------------------------------------------------

###-----Building Stack
cov_bld <- datt[,c("x3","x15","x27","x30", "x31","x33","x34", "x35","x36",
                   "x38","x46","x50","x52","set_prov", "set_typ", 
                   "prov", "IDsp")]; dim(cov_bld)

#---Build the stack for the bld
datt$BLDG21 <- datt$BLDG21 + 1# transformed sto have at least a value of 1 
stk_bld <- inla.stack(data=list(y=log(datt$BLDG21)), #the response
                      
                      A=list(A,1),  #the A matrix; the 1 is included to make the list(covariates)
                      
                      effects=list(c(list(Intercept=1), #the Intercept
                                     iset),  #the spatial index
                                   #the covariates
                                   list(cov_bld)
                      ), 
                      #this is a quick name so you can call upon easily
                      tag='est_bld')


##-----------------------------------------------------
# Model fitting
#-------------------------------------------------------

# Model for building intensity (this is the best fit model following initial multiple iterations)
fbld3<- y ~ -1 + Intercept +  x3 + x15 + x27 + x30 + x31 + x33 + x34 + x35 + x36 + 
  x38 + x46 + x50 + x52 + f(spatial.field, model=spde) + f(IDsp, model='iid') + f(set_prov, model="iid")

bmod3<-inla(fbld3, #the formula
            data=inla.stack.data(stk_bld,spde=spde),  #the data stack
            family= 'gaussian',   #which family the data comes from
            control.predictor=list(A=inla.stack.A(stk_bld),compute=TRUE),  #compute gives you the marginals of the linear predictor
            control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
            verbose = FALSE) #can include verbose=TRUE to see the log of the model runs
summary(bmod3) #----model summary
bmod3$summary.fix #--extract fixed effects 
bmod3$summary.random #--extract random effects 
bind3 <-inla.stack.index(stk_bld, "est_bld")$data #--estimation indices 
bfit3 <- exp(bmod3$summary.linear.predictor[bind3,"mean"]) #--extract the backtransformed building intensity
sum(bfit3)
bfit3L <- exp(bmod3$summary.linear.predictor[bind3,"0.025quant"]) #--extract the backtransformed lower building intensity
bfit3U <- exp(bmod3$summary.linear.predictor[bind3,"0.975quant"]) #--extract the backtransformed upper building intensity


# explore model outputs
(bld3<- cbind(datt$BLDG21, bfit3))
apply(bld3, 2, sum, na.rm=T)
plot(datt$BLDG21, bfit3, col=c(1,2))
abline(a=0, b=1)
cor(datt$BLDG21, bfit3)


# Add the predicted building 
datt$bld_pred <- bfit3

#--extract fixed effects 
betab <- bmod3$summary.fix 


##---Carry out model checks using waic
mod.fitb<- bmod3$waic$waic
(b_mets <- mod_metrics2(datt$BLDG21, bfit3, bmod3))


#------------------------------------------------
#  Define density
#------------------------------------------------
# TSBHM - uses the bias-adjusted settlement/building data
datt$dens2 <- datt$POPN/datt$bld_pred 
hist(log(datt$dens2), col="brown")
datt$dens2[datt$dens2==0] = 0.000001 ####4918 16520 16521 16527 18142 18149 18262 18267 18297 18414 20900 20908; 12 CU's no observations


# BHM  - uses the imperfectly observed settlement/building data directly 
datt$bld <- datt$BLDG21 # rename the building intensity variable
datt$bld[datt$POPN==0] =NA # set all 0 POPULATION to NA to allow for prediction and avoid dividing by zero
datt$dens3 <- datt$POPN/datt$bld
datt$dens3[datt$dens3==0] = 0.000001 ####4918 16520 16521 16527 18142 18149 18262 18267 18297 18414 20900 20908; 12 CU's no observations

####---Density Stack
covars_dens <- datt[,c("x3","x11","x15","x29","x30", "x31","x32","x35", "x36",
                       "x39","x45","x46","x48","x51","x52","set_prov", "set_typ", 
                       "prov", "IDsp")]; dim(covars_dens)

#---Build the stack for the training set
##---Bias Corrected (TSBHM)
stk_dens <- inla.stack(data=list(y=datt$dens2), #the response
                       
                       A=list(A,1),  #the A matrix; the 1 is included to make the list(covariates)
                       
                       effects=list(c(list(Intercept=1), #the Intercept
                                      iset),  #the spatial index
                                    #the covariates
                                    list(covars_dens)
                       ), 
                       #this is a quick name so you can call upon easily
                       tag='est_dens')


##--(BHM)
stk3<- inla.stack(data=list(y=datt$dens3), #the response
                  
                  A=list(A,1),  #the A matrix; the 1 is included to make the list(covariates)
                  
                  effects=list(c(list(Intercept=1), #the Intercept
                                 iset),  #the spatial index
                               #the covariates
                               list(covars_dens)
                  ), 
                  #this is a quick name so you can call upon easily
                  tag='est3')


###---SPATIAL---------------------------------
##========================== TSBHM
form1a <- y ~ -1 + Intercept +  x3 + x11 + x15 + x29 + x30 + x31 + x32 + x35 + x36 + 
  x39 + x45+ x46 + x48 + x51 + x52 + f(spatial.field, model=spde) + f(IDsp, model='iid') +
  f(set_typ, model='iid')

mod1a <-inla(form1a, #the formula
             data=inla.stack.data(stk_dens,spde=spde),  #the data stack
             family= 'gamma',   #which family the data comes from
             control.predictor=list(A=inla.stack.A(stk_dens),compute=TRUE),  #compute gives you the marginals of the linear predictor
             control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
             verbose = FALSE) #can include verbose=TRUE to see the log of the model runs
#summary(mod1c) #----model summary


betagg <- mod1a$summary.fix #--extract fixed effects 
#write.csv(round(betagg,5), paste0(results_path,"/updated2/Betas_for_density_model.csv"))

#mod1$summary.random #--extract random effects 
ind1 <-inla.stack.index(stk_dens, "est_dens")$data #--estimation indices 

# predicted population density and uncertainties
fit1a <- exp(mod1a$summary.linear.predictor[ind1,"mean"]) #--extract the backtransformed pop_hat
fit1aU <- exp(mod1a$summary.linear.predictor[ind1,"0.975quant"]) #--extract the backtransformed pop_hat
fit1aL <- exp(mod1a$summary.linear.predictor[ind1,"0.025quant"]) #--extract the backtransformed pop_hat


# predicted population count and uncertainties
fit11a <- fit1a*datt$bld_pred
fit11aU <- fit1aU*datt$bld_pred
fit11aL <- fit1aL*datt$bld_pred
(sum_fit <- sum(fit1a*datt$bld_pred))
(POPa <- cbind(datt$POPN, fit11a))
apply(POPa, 2, sum, na.rm=T)
plot(datt$POPN, fit11a, col=c("red", "blue"))
abline(a=0, b=1)
cor(datt$POPN, fit11a)

#######
###---SPATIAL---------------------------------
##========================== BHM
form1b <- y ~ -1 + Intercept +  x3 + x11 + x15 + x29 + x30 + x31 + x32 + x35 + x36 + 
  x39 + x45+ x46 + x48 + x51 + x52 + f(spatial.field, model=spde) + f(IDsp, model='iid') +
  f(set_typ, model='iid')

mod1b <-inla(form1b, #the formula
             data=inla.stack.data(stk3,spde=spde),  #the data stack
             family= 'gamma',   #which family the data comes from
             control.predictor=list(A=inla.stack.A(stk3),compute=TRUE),  #compute gives you the marginals of the linear predictor
             control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
             verbose = FALSE) #can include verbose=TRUE to see the log of the model runs
#summary(mod1c) #----model summary

#mod1$summary.random #--extract random effects 
ind3 <-inla.stack.index(stk3, "est3")$data #--estimation indices 
fit1b <- exp(mod1b$summary.linear.predictor[ind3,"mean"]) #--extract the backtransformed pop_hat
fit1bU <- exp(mod1b$summary.linear.predictor[ind3,"0.975quant"]) #--extract the backtransformed pop_hat
fit1bL <- exp(mod1b$summary.linear.predictor[ind3,"0.025quant"]) #--extract the backtransformed pop_hat

fit11b <- fit1b*datt$bld
fit11bU <- fit1bU*datt$bld
fit11bL <- fit1bL*datt$bld
(sum_fit <- sum(fit1b*datt$bld, na.rm=T))
(POPb <- cbind(datt$POPN, fit11b))
apply(POPb, 2, sum, na.rm=T)
plot(datt$POPN, fit11b, col=c("red", "blue"))
abline(a=0, b=1)
cor(datt$POPN, fit11b)



(met1a <- mod_metrics2(datt$POPN, fit11a, mod1a))
(met1b <- mod_metrics2(datt$POPN, fit11b, mod1b))


(met.1 <- data.frame(mod1a = unlist(met1a),
                     mod1b = unlist(met1b)))


par(mfrow=c(1,2))
plot(datt$POPN, fit11a, col=c("red", "blue"))
abline(a=0, b=1)
plot(datt$POPN, fit11b, col=c("red", "blue"))
abline(a=0, b=1)
par(mfrow=c(1,1))




#Compute statistics in terms or range and variance
post_effa <- inla.spde2.result(inla = mod1a, name = "spatial.field",
                               spde = spde, do.transf = TRUE)

post_effb <- inla.spde2.result(inla = mod1b, name = "spatial.field",
                               spde = spde, do.transf = TRUE)

out_path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP000008_UNFPA_PNG/Working/Chris/paper"
##--TSBHM----Posterior parameter estimates
beta_a <- round(mod1a$summary.fix[,c(1,2,3,5)],4) #--extract fixed effects 
prec_gam_a <- round(mod1a$summary.hyperpar[c(1,4,5),c(1,2,3,5)],4)#--Precision for Gamma observations AND RANDOM EFFECTS
marg_var_a <- round(as.vector(unlist(inla.zmarginal(post_effa$marginals.variance.nominal[[1]]))[c(1,2,3,5)]),4)
waic_a <- mod1a$waic$waic
cpo_a <- -sum(log(mod1a$cpo$cpo), na.rm=T)
post_est_a <- bind_rows(beta_a, prec_gam_a)
rownames(post_est_a) <- c("Int",rownames(post_est_a)[2:16], "eps", "tau_iid", "tau_st")
write.csv(post_est_a, paste0(out_path, "/post_est_mod1a.csv"))



##--BHM----Posterior parameter estimates
beta_b <- round(mod1b$summary.fix[,c(1,2,3,5)],4) #--extract fixed effects 
prec_gam_b <- round(mod1b$summary.hyperpar[c(1,4,5),c(1,2,3,5)],4)#--Precision for Gamma observations AND RANDOM EFFECTS
marg_var_b <- round(as.vector(unlist(inla.zmarginal(post_effb$marginals.variance.nominal[[1]]))[c(1,2,3,5)]),4)
waic_b <- mod1b$waic$waic
cpo_b <- -sum(log(mod1b$cpo$cpo), na.rm=T)
post_est_b <- bind_rows(beta_b, prec_gam_b)
rownames(post_est_b) <- c("Int",rownames(post_est_b)[2:16], "eps", "tau_iid", "tau_st")
write.csv(post_est_b, paste0(out_path, "/post_est_mod1b.csv"))


var_b <- inla.zmarginal(post_effb$marginals.variance.nominal[[1]])



write.csv(round(betab,5), paste0(results_path,"/updated2/Betas_for_building_intensity_model.csv"))


#summary(mod1c) #----model summary
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#####-----POSTERIOR SIMULATION-----------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

########################

stypp <- function(mod, dat)
{
  #datt$set_typ2 <- rep(1, nrow(datt))
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
datt$set_typ_a <- stypp(mod1a, datt)
datt$set_typ_b <- stypp(mod1b, datt)

datt4 <- datt
datt$set_re <- datt$set_typ_a 
datt4$set_re <- datt$set_typ_b 
datt$BLD <- datt$bld_pred
datt4$BLD <- datt$bld
names(datt4)

###################################=======-------------------------GG
simPops_gg <- function(model, dat, Aprediction, run)
{
  fixedeff  <- dens_hat <- pop_hat <- matrix(0, nrow=nrow(dat), ncol = run)
  # inla.seed = as.integer(runif(1)*.Machine$integer.max)
  inla.seed =  481561959
  set.seed(inla.seed)
  print(inla.seed)
  m1.samp <- inla.posterior.sample(run, model, seed = inla.seed ,selection=list(x3=1, x11=1, x15=1,
                                                                                x29=1, x30=1, x31=1,
                                                                                x32=1, x35=1, x36=1,
                                                                                x39=1, x45=1, x46=1,
                                                                                x48=1, x51=1, x52=1),num.threads="1:1")
  
  sfield_nodes_mean <- model$summary.random$spatial.field['mean']
  field_mean <- (Aprediction%*% as.data.frame(sfield_nodes_mean)[, 1])
  for(i in 1:run)
  {
    #fixedeff[,i] <- model$summary.fixed['Intercept', 'mean'] +
    fixedeff[,i] <-  
      #m1.samp[[i]]$latent[1,] +
      model$summary.fixed['Intercept', 'mean'] +
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
      model$summary.random$IDsp['mean'][,1] +
      dat$set_re +
      #rnorm(nrow(dat), 0, 1/m1.samp[[i]]$hyperpar[4]) + 
      field_mean[,1]
    
    dens_hat[,i]<- exp(fixedeff[,i])
    pop_hat[,i] <- dens_hat[,i]*dat$BLD
  }
  
  #mean_pop_hat1 <- dat$pop_hat1 #
  mean_dens_hat <- apply(dens_hat, 1, mean, na.rm=T) #
  mean_pop_hat <- apply(pop_hat, 1, mean, na.rm=T) #
  median_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.5), na.rm=T) #
  lower_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.025), na.rm=T) #
  upper_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.975), na.rm=T) #
  sd_pop_hat <- apply(pop_hat, 1, sd, na.rm=T) #
  uncert_pop_hat <- (upper_pop_hat - lower_pop_hat)/mean_pop_hat#
  
  
  dat$mean_dens_hat <- mean_dens_hat
  dat$mean_pop_hat <- mean_pop_hat
  dat$median_pop_hat <- median_pop_hat
  dat$lower_pop_hat <- lower_pop_hat
  dat$upper_pop_hat <- upper_pop_hat
  dat$uncert_pop_hat <- uncert_pop_hat
  dat$sd_pop_hat <- sd_pop_hat
  
  
  
  output <- list(pop_hat = pop_hat,
                 est_data = dat)
  
}

#rm(dat.sim, sim.pops_nb, sim.pops, llg.est, prov.est2)
run=2000
#run=100
system.time(str(sim.pops_gg <- simPops_gg(mod1a,datt, A, run)))  # TSBHM
system.time(str(sim.pops_gg2 <- simPops_gg(mod1b,datt4, A, run))) # BHM


# Visualise the posterior samples 
# Load ggplot2
library(ggplot2)

##-----Traceplots 
samp <- sample(nrow(datt), 6)
par(mfrow=c(3,2), mar=c(5,5,2,1))
for(j in samp)
{
  plot.ts(sim.pops_gg$pop_hat[j,], ylab = "pop_hat", col="steel blue",
          cex.main=2,  main=paste0("Pop_hat samples for CU ", datt$CU_Name[j], sep=""))
  abline(h=mean(sim.pops_gg$pop_hat[j,]), lwd=2, col=2)
}



#####----Join the posterior draws to the dataframe------------------------------------------------------
dim(datt)
names(datt)
data.all <- datt[,c("Prov_Name" ,"Dist_Name", "LLG_Name", "Ward_Name")]

# TSBHM
data.all1a <- data.frame(cbind(data.all, sim.pops_gg$pop_hat))#---Contains all the simulated posterior matrix
dim(data.all1a);names(data.all1a)


# BHM
data.all1b<- data.frame(cbind(data.all, sim.pops_gg2$pop_hat))#---Contains all the simulated posterior matrix
dim(data.all1b);names(data.all1b)


####-----------ADMIN TOTALS AND UNCERTAINTIES----------------------------------------
thin=1 #----Thinning choses every 5th sample for posterior inference 
thinn <- 4+(1:run)#--There is a burnin period of the first 20% of the total sample


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
  
  return(estimates <- data.frame(estimates=unlist(list(total=tot_mean, lower=tot_lower, median=tot_median, upper=tot_upper))))
}
(national1 <- nat_total(data.all1a, thinn)) # TSBHM
(national2<- nat_total(data.all1b, thinn)) # BHM



write.csv(national1, file=paste0(results_path, "/updated2/National_estimates_main_gamma-gaussian.csv"))


###########
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
#(prov.est_gg <- prov_est_gg(data.all,thinn))

prov1 <- prov_est_gg(data.all1a,thinn) # TSBHM
prov2<- prov_est_gg(data.all1b,thinn) # BHM
write.csv(prov1, file=paste0(res_path, "/Province_data_TSBHM.csv"), row.names=F)
write.csv(prov2, file=paste0(res_path, "/Province_data_BHM.csv"), row.names=F)


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

(dist1 <- dist_est_gg(data.all1a,thinn)) # TSBHM
(dist2<- dist_est_gg(data.all1b,thinn)) # BHM
write.csv(dist1, file=paste0(res_path, "/District_data_TSBHM.csv"), row.names=F)
write.csv(dist2, file=paste0(res_path, "/District_data_BHM.csv"), row.names=F)


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

llg1 <- llg_est_gg(data.all1a,  thinn) # TSBHM 
llg2<- llg_est_gg(data.all1b,  thinn) # BHM 


write.csv(llg1, file=paste0(res_path, "/LLG_data_TSBHM.csv"), row.names=F)
write.csv(llg2, file=paste0(res_path, "/LLG_data_BHM.csv"), row.names=F)


### --- Ward level estimates with uncertainties 
wd_est_gg <- function(datl, thinn)
{
  lnames <- unique(datl$Ward_Name)
  outl <- matrix(0, nrow=length(lnames), ncol=3)
  
  for(j in 1:length(lnames))
  {
    wd <- datl[datl$Ward_Name==lnames[j],]
    ltots <- apply(wd[,thinn], 2, sum, na.rm=T)
    ltots_sd <- sd(ltots, na.rm=T)
    
    ltot_mean1  <- mean(ltots, na.rm=T)
    ltot_lower <- quantile(ltots, probs=c(0.025))
    ltot_upper <- quantile(ltots, probs=c(0.975))
    
    lestimates <- round(c(ltot_mean1, ltot_lower, ltot_upper), 3)
    outl[j,] <- lestimates
  }
  outl <- data.frame(outl)
  return(wd_est <- data.frame(names = lnames,
                               total = outl[,1],
                               lower = outl[,2],
                               upper = outl[,3]))
}

wd1 <- wd_est_gg(data.all1a,  thinn) # TSBHM 
wd2<- wd_est_gg(data.all1b,  thinn) # BHM 


write.csv(wd1, file=paste0(res_path, "/Ward_data_TSBHM.csv"), row.names=F)
write.csv(wd2, file=paste0(res_path, "/Ward_data_BHM.csv"), row.names=F)

#------------------------------------------------------------
###---CU level estimates
#----------------------------------------------------------
cu1 <- sim.pops_gg$est_data # TSBHM 
cu2 <- sim.pops_gg2$est_data # BHM 

write.csv(cu1, file=paste0(res_path, "/cu_data_TSBHM.csv"), row.names=F)
write.csv(cu2, file=paste0(res_path, "/cu_data_BHM.csv"), row.names=F)


#------------------------------------------
# Extract psoterior data for each appraoch
var2Add <- c("Prov_Name", "Dist_Name", "LLG_Name", "Ward_Name", "CU_Name", "POPN", "BLDG21", "bld_pred",
             "set_typ", "mean_dens_hat", "mean_pop_hat", "lower_pop_hat", "upper_pop_hat", "uncert_pop_hat")

# TSBHM
dim(dat_tsbhm <- cu1) 
dat_tsbhm$method <- rep("TSBHM", nrow(dat_tsbhm))
dat_tsbhm$ldens <- log(dat_tsbhm$mean_dens_hat) # log posterior pop density


# BHM
dim(dat_bhm <-cu2) 
dat_bhm$method <- rep("BHM", nrow(dat_bhm))
dat_bhm$ldens <- log(dat_bhm$mean_dens_hat)  # log posterior pop density


# Combined data

dim(dat_all <- rbind(dat_tsbhm, dat_bhm))

#--density plot
pd <- ggdensity(dat_all, x = "ldens",
                add = "mean", rug = TRUE,
                color = "method", fill = "method"#,
                #palette = c("#00AFBB", "#E7B800")
) 
pd

pdens <- ggpar(pd, xlab="Log of Predicted Population Density", ylab="Density",
               legend = "top", legend.title = "Method",size=22,
               font.legend=c(18),
               palette = "jco",
               font.label = list(size = 15, face = "bold", color ="red"),
               font.x = c(16),
               font.y = c(16),
               xtickslab.rt = 45, 
               ytickslab.rt = 45)
pdens

##---Histogram plots

gh <- gghistogram(dat_all, x = "ldens",
                  add = "mean", rug = TRUE,
                  size=2,
                  color = "method", fill = "method", bins=50#,
                  # palette = c("#00AFBB", "#E7B800")
) 

phdens <- ggpar(gh, xlab="Log of Predicted Population Density", ylab="Count",
                legend = "top", legend.title = "Method",size=22,
                font.legend=c(20),
                palette = c("jco"),
                #[1] "#FF0000" "#FF1C00" "#FF3800" "#FF5500" "#FF7100" "#FF8D00" "#FFAA00", "#E5F5F9","#99D8C9" 
                #[8] "#FFC600" "#FFE200" "#FFFF00"
                #"#E5F5F9" "#99D8C9" "#FFE200"
                #colour = "bld_cover",
                #shape= "bld_cover",
                yscale = c("none"),
                font.label = list(size = 18, face = "bold", color ="red"),
                font.x = c(22),
                font.y = c(20),
                font.main=c(14),
                font.xtickslab =c(18),
                font.ytickslab =c(20),
               # ylim= c(0, max(p100m$pop)),
                xtickslab.rt = 45, ytickslab.rt = 45)
phdens


###Scatter plots: CU-level
names(dat_all)
dat_all$method <- factor(dat_all$method)
pphat <- dat_all %>%
  ggplot(aes(x=POPN, y=mean_pop_hat))+
  # geom_point(aes(colour=method))+
  geom_pointrange(aes(ymin = lower_pop_hat, ymax = upper_pop_hat, colour=method), size=0.5)+
  #geom_errorbar(aes(ymin = lower_pop_hat, ymax = upper_pop_hat, colour=method), size=1,width = 0.5)+
  geom_smooth(aes(colour=method),method="lm", se=F)+
  theme_bw()+
  theme(strip.text = element_text(size = 15),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14),
        #panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  facet_wrap(~method)
pphat
pv <- ggpar(pphat, xlab="Observed population(count)", ylab="Predicted population(count)",
            legend = "top", legend.title = "Method",
            palette = "jco",
            font.label = list(size = 15, face = "bold", color ="red"),
            font.x = c(16),
            font.y = c(16),
            xtickslab.rt = 45, 
            ylim= c(0, max(dat_all$POPN, na.rm=T)),
            ytickslab.rt = 45)

pv




# install.packages("remotes")
#remotes::install_github("R-CoderDotCom/ridgeline@main")

#library(ridgeline)

#ridgeline(chickwts$weight, chickwts$feed)

#library(ridgeline)

##

## Calculate  Error rates 
# TSBHM
dim(cu1b <- cu1 %>% drop_na(POPN)) 
predTot1 <- sum(cu1b$mean_pop_hat, na.rm=T)
trueTot1 <- sum(cu1b$POPN,na.rm=T)
#bb1 <- sum(cu1b$mean_pop_hat-cu1b$POPN,na.rm=T)/nrow(cu1)
bias1b <- abs(predTot1-trueTot1)/trueTot1


# BHM
dim(cu2b <- cu2 %>% drop_na(POPN)) 
predTot2 <- sum(cu2b$mean_pop_hat, na.rm=T)
trueTot2 <- sum(cu2b$POPN,na.rm=T)
#bb2 <- sum(cu2b$mean_pop_hat-cu2b$POPN,na.rm=T)/nrow(cu2)
bias2b <- abs(predTot2-trueTot2)/trueTot2

## Calculate  relative error rates
(rbiasb <- (bias1b/bias2b))  
#(rbiasb2 <- (bb1/bb2)) 
## Calculate percentage reduction in error rates
(1-rbiasb)*100
hist(log(cu2$POPN-cu2$mean_dens_hat))



require(ggplot2)
ggplot(cu1, aes(x="mean_pop_hat"))+
geom_boxplot()

sum(cu.est$mean, na.rm=T)
Vars2Include <- c("CLUSTER_ID", "Uniq_ID", "Prov_Name" ,"Dist_Name", "LLG_Name", "Ward_Name",
                  "CU_Name", "TYPE", "lon", "lat", "SOURCE", "POPN", "BLDG21", "mean", "median",
                  "lower", "upper", "uncertainty", "mean_dens_hat")
names(cu1)
dim(cu.data1 <- cu1[,Vars2Include]); head(cu1)
cu.data1$method <- rep("TSBHM", nrow(cu.data1))
dim(cu.data2 <- cu2[,Vars2Include]); head(cu2)
cu.data2$method <- rep("BHM", nrow(cu.data2))


dim(cu_dat <- rbind(cu.data1,cu.data2))

## ward summary
(wd.data1 <- cu.data1 %>% group_by(Ward_Name) %>%
  summarise(mean = sum(mean, na.rm=T),
            pop = sum(POPN, na.rm=T)))
wd.data1$method <- rep("TSBHM", nrow(wd.data1))


## lga summary
(lg.data1 <- cu.data1 %>% group_by(LLG_Name) %>%
    summarise(mean = sum(mean, na.rm=T),
              pop = sum(POPN, na.rm=T)))
lg.data1$method <- rep("TSBHM", nrow(lg.data1))



##  district summary
(dst.data1 <- cu.data1 %>% group_by(Dist_Name) %>%
    summarise(mean = sum(mean, na.rm=T),
              pop = sum(POPN, na.rm=T)))
dst.data1$method <- rep("TSBHM", nrow(dst.data1))



##  province summary
(prv.data1 <- cu.data1 %>% group_by(Prov_Name) %>%
    summarise(mean = sum(mean, na.rm=T),
              pop = sum(POPN, na.rm=T)))
prv.data1$method <- rep("TSBHM", nrow(prv.data1))





####
##
ppt <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP000008_UNFPA_PNG/Working/Amy/PNG_report/FINAL_models_310123/Shapefiles/gamma_gaussian"
##
shpp <- st_read(paste0(ppt, "/gamma_gaussian_CU.shp"))
names(shpp)
plot(shpp["AREA_KM2"])
names(shp.cu)

plot(shp.cu["obs"])

library(tmap)
library(tmaptools)

tmap_options(check.and.fix = TRUE)
tmap_mode(mode = "plot")

#require(rgdal)
crs.UTM = CRS("+proj=utm +zone=55 +datum=WGS84 +units=m +no_defs")#png is zone 55
shpcu_UTM2 = st_transform(shp.cu, crs.UTM)
#--Observed counts


##--predicted based on best fit gg model
pred1 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="pop_hat1", title="Predicted Count",
              legend.hist=F, 
              palette=viridis(100), 
              legend.show =T,
              breaks=c(1,100,300,600, 1000,3000, 7000, 12000)
              )+
  tm_borders(col="white")+
  tm_layout(legend.outside = F, legend.text.size=1.5, legend.title.size=2)+
  tm_compass(position = c("right", "bottom"), text.size=1.5)+
  tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "", frame=FALSE) 



pred2 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="pred2", title="Predicted Count",
              legend.hist=F, 
              palette="plasma", 
              legend.show =F,
              breaks=c(1,100,300,600, 1000,3000, 7000, 12000)
  )+
  tm_layout(legend.outside = F, legend.text.size=1.5, legend.title.size=2)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "", frame=FALSE) 


##--predicted based on best fit gg model
uncert1 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="Uncertaint", title="Uncertainty",
              legend.hist=F, 
              palette=viridis(100), 
              legend.show = T,
              breaks=c(0,0.1,0.2, 0.3,0.4,0.5,0.6,1.2,2.5,4))+
  tm_borders(col="white")+
  tm_layout(legend.outside = F, legend.text.size=1.5, legend.title.size=2)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1.5, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "", frame=FALSE) 




class(cu_dat$method <- factor(cu_dat$method))
table(cu_dat$method)

# Scatter plots (sp) with regression linem
library(ggpubr)
pc <- ggscatter(cu_dat, x = "POPN", y = "mean",
                add = "reg.line",               # Add regression line
                conf.int = TRUE,                # Add confidence interval
                color = "method", palette = "jco", # Color by groups "cyl"
                #shape = "Method",                   # Change point shape by groups "cyl"
                facet.by ="method"
)

pcb <- ggpar(pc, xlab="Observed Population", ylab="Predicted Population",
             legend = "right", legend.title = "Area Categories",
             font.label = list(size = 15, face = "bold", color ="red, blue"),
             font.x = c(16),
             font.y = c(16),
             #xlim=c(0,10000),
             xtickslab.rt = 45, ytickslab.rt = 45)

pcb




# CU level scatter plots 

pc <- ggplot(cu_dat, aes(POPN, mean, colour=method)) +
#geom_pointrange(aes(ymin = lower, ymax = upper))+
  #geom_smooth(method="lm")+
  geom_linerange(aes(ymin = lower, ymax = upper), size=1) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size=0.5)+
  geom_errorbar(aes(ymin = lower, ymax = upper), size=1,width = 0.5)+
  geom_line(aes(group = method), size=1) +
 theme_bw()+
  #facet_wrap(~method, scales="free", nrow=1)+
  theme(strip.text = element_text(size = 15),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14))#+
  #scale_fill_manual(values=c("#00AFBB", "#E7B800"))
  

pcb <- ggpar(pc, xlab="Observed Population", ylab="Predicted Population",
             legend = "right", legend.title=element_text("Method",size=22),
             font.legend=c(18),
             palette = c("#0D0887FF", "#00AFBB"),
             font.label = list(size = 18, face = "bold", color ="red"),
             font.x = c(22),
             font.y = c(18),
             font.main=c(14),
             font.xtickslab =c(18),
             font.ytickslab =c(16),
             ylim=c(0, 13000),
             xtickslab.rt = 45, ytickslab.rt = 45)

pcb


#####
cu_dat$logdens <- log(cu_dat$mean_dens_hat)
gh <- gghistogram(cu_dat, x = "logdens",
                  add = "mean", rug = TRUE,
                  color = "method", fill = "method", bins=50) 

ph <- ggpar(gh, xlab="Log of Predicted Population Density", ylab="Count",
            legend = "right", legend.title=element_text("Method",size=22),
            font.legend=c(18),
            palette = c("#0D0887FF", "#00AFBB"),
            font.label = list(size = 18, face = "bold", color ="red"),
            font.x = c(22),
            font.y = c(18),
            font.main=c(14),
            font.xtickslab =c(18),
            font.ytickslab =c(16),
            #ylim=c(0, 13000),
            xtickslab.rt = 45, ytickslab.rt = 45)
ph


# CU level scatter plots 

pc <- ggplot(cu_dat, aes(POPN, mean, colour=method)) +
  #geom_pointrange(aes(ymin = lower, ymax = upper))+
  #geom_smooth(method="lm")+
  geom_linerange(aes(ymin = lower, ymax = upper), size=1) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size=0.5)+
  geom_errorbar(aes(ymin = lower, ymax = upper), size=1,width = 0.5)+
  geom_line(aes(group = method), size=1) +
  theme_bw()+
  #facet_wrap(~method, scales="free", nrow=1)+
  theme(strip.text = element_text(size = 15),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14))#+





#######

####

# Basic plot
# +++++++++++++++++++++++++++

cu <- ggscatter(cu.data1, x = "POPN", y = "mean",
                color = "black", shape = 16, size = 2, # Points color, shape and size
                add = "reg.line",  # Add regressin line
                #ellipse = TRUE,
                add.params = list(color = "red", fill = "light green"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)

cu11 <- ggpar(cu, xlab="Observed population count", ylab="Predicted population count",
              legend = "right", legend.title=element_text("Method",size=22),
              font.legend=c(18),
              palette = c("#0D0887FF"),
              font.label = list(size = 18, face = "bold", color ="red"),
              font.x = c(22),
              font.y = c(18),
              font.main=c(14),
              font.xtickslab =c(18),
              font.ytickslab =c(16),
              #ylim=c(0, 13000),
              xtickslab.rt = 45, ytickslab.rt = 45)
cu11




wd <- ggscatter(wd.data1, x = "pop", y = "total",
                color = "black", shape = 16, size = 2, # Points color, shape and size
                add = "reg.line",  # Add regressin line
               # ellipse = TRUE,
                add.params = list(color = "red", fill = "light green"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)

wd11 <- ggpar(wd, xlab="Observed population count", ylab="Predicted population count",
              legend = "right", legend.title=element_text("Method",size=22),
              font.legend=c(18),
              palette = c("00AFBB"),
              font.label = list(size = 18, face = "bold", color ="red"),
              font.x = c(22),
              font.y = c(18),
              font.main=c(14),
              font.xtickslab =c(18),
              font.ytickslab =c(16),
              #ylim=c(0, 13000),
              xtickslab.rt = 45, ytickslab.rt = 45)
wd11



#FF3800"

lg <- ggscatter(lg.data1, x = "pop", y = "mean",
                color = "black", shape = 16, size = 2, # Points color, shape and size
                add = "reg.line",  # Add regressin line
                #ellipse = TRUE,
                add.params = list(color = "red", fill = "light green"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)
lg11 <- ggpar(lg, xlab="Observed population count", ylab="Predicted population count",
              legend = "right", legend.title=element_text("Method",size=22),
              font.legend=c(18),
              palette = c("black"),
              font.label = list(size = 18, face = "bold", color ="red"),
              font.x = c(22),
              font.y = c(18),
              font.main=c(14),
              font.xtickslab =c(18),
              font.ytickslab =c(16),
              #ylim=c(0, 13000),
              xtickslab.rt = 45, ytickslab.rt = 45)
lg11




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


ggarrange(cu11, wd11, lg11, dst11, prv11, nrow=3, ncol=2)

# install.packages("MASS")
library(MASS)


cu.data1b <- cu.data1 %>% drop_na(POPN, mean)
kern <- kde2d(cu.data1$POPN[!is.na(cu.data1$POPN)], cu.data1$mean[!is.na(cu.data1$POPN)])
contour(kern, drawlabels = FALSE, nlevels = 6,
        col = rev(heat.colors(6)), add = TRUE, lwd = 3)




#---------------------------------------------------------------------------
library(ggplot2)

commonTheme = list(labs(color="Density",fill="Density",
                        x="RNA-seq Expression",
                        y="Microarray Expression"),
                   theme_bw(),
                   theme(legend.position=c(0,1),
                         legend.justification=c(0,1)))

names(prov1)
prvd1 <- prov1 %>% select(c("total", "lower")) %>% drop_na()
cud1 <- cud1[,c("Predicted", "Observed")]
cud2 <- cu_sc_lng %>% filter(Method=="TSBHM") %>% drop_na()


names(cu1)

# Dist_Name
# Prov_Name
# Ward_Name 
# CU_Name
# LLG_Name
lgdt <- cu1 %>% group_by(Dist_Name)%>% drop_na(mean) %>%
  summarise(Obs = sum(POPN, na.rm=T),
            Pred = round(sum(mean, na.rm=T)))


lgdt <- cu2
lgdt1 <- lgdt[,c("POPN", "mean")]
#lgdt1 <- lgdt[,c("Obs", "Pred")]
df1 = data.frame(lgdt1); colnames(df1) = c("x","y")

commonTheme = list(labs(color="Density",fill="Density",
                        x="Observed Population",
                        y="Predicted Popualtion"),
                   theme_bw(),
                   theme(legend.position=c(0,1),
                         legend.justification=c(0,1)))

dens1 <- ggplot(data=lg.data1,aes(pop,mean)) + 
 stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='blue') + 
 # stat_density2d(geom='polygon',colour='blue') + 
  #scale_fill_continuous(low="green",high="red") +
  scale_fill_continuous(type = "viridis") +
  geom_smooth(method=lm,linetype=2,colour="red",se=F) + 
  guides(alpha="none") +
  #scale_y_continuous(breaks=seq(0,max(df1$y, na.rm=T),50000))+
 # scale_x_continuous(breaks=seq(0,max(df1$x, na.rm=T),50000))+
  
  
 # scale_y_continuous(breaks=seq(0,max(df1$y, na.rm=T),50000))+
  #scale_x_continuous(breaks=seq(0,max(df1$x, na.rm=T),50000))+
  geom_point() + commonTheme


####

plot_dens2<- ggpar(dens1, ylab="Predicted Population", xlab="Observed Population",
                 legend = "right", legend.title=element_text("Density",size=22) ,
                 font.legend=c(20),
                 legend.show=F,
                 # legend.text=element_text(size=22),
                 #main ="Lollipop plot of percentage relative  bias reduction by TSBHM",
                 #palette = c("#00AFBB","blue", "red"),
                 #[1] "#FF0000" "#FF1C00" "#FF3800" "#FF5500" "#FF7100" "#FF8D00" "#FFAA00", "#E5F5F9","#99D8C9" 
                 #[8] "#FFC600" "#FFE200" "#FFFF00"
                 #"#E5F5F9" "#99D8C9" "#FFE200"
                 #colour = "bld_cover",
                 #shape= "bld_cover",
                 #scale_x_continuous(limits=c(0,100), breaks=c(20, 40,60, 80, 100)),
                 font.label = list(size = 18, face = "bold", color ="red"),
                 font.x = c(20),
                 font.y = c(20),
                 font.main=c(14),
                 font.xtickslab =c(18),
                 font.ytickslab =c(18),
                 #xticks.by = TRUE,
                 #ylim= c(0, max(p100m$pop)),
                 xtickslab.rt = 80, ytickslab.rt = 45)

plot_dens2



# ---------------------------------------------------------------------------------
## CU level density and scatter plots
# Bin size control + color palette
lgdt <- cu1 # TSBHM
lgdt1 <- lgdt[,c("POPN", "mean")]
#lgdt1 <- lgdt[,c("Obs", "Pred")]
df1 = data.frame(lgdt1); colnames(df1) = c("x","y")

cor(df1$x[!is.na(df1$x)], df1$y[!is.na(df1$x)])
cup1 <- ggplot(df1, aes(x=x, y=y) ) + 
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  #geom_point()+
  geom_smooth(method=lm,linetype=2,colour="red",se=F) + 
  theme_bw()

plot_cup1<- ggpar(cup1, ylab="Predicted Population", xlab="Observed Population",
                  legend = "right", legend.title=element_text("Count",size=22) ,
                  font.legend=c(20),
                  legend.show=F,
                  font.label = list(size = 20, face = "bold", color ="red"),
                  font.x = c(20),
                  font.y = c(20),
                  font.main=c(16),
                  font.xtickslab =c(20),
                  font.ytickslab =c(20),
                  #xticks.by = TRUE,
                  xlim= c(0, max(df1$x)),
                   xtickslab.rt = 45, ytickslab.rt = 45)

plot_cup1


#### - BHM 
lgdt2 <- cu2
lgdt2 <- lgdt2[,c("POPN", "mean")]
#lgdt1 <- lgdt[,c("Obs", "Pred")]
df2 = data.frame(lgdt2); colnames(df2) = c("x","y")

cup2 <- ggplot(df2, aes(x=x, y=y) ) +
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  geom_point()+
  geom_smooth(method=lm,linetype=2,colour="red",se=F) + 
  theme_bw()

plot_cup2<- ggpar(cup2, ylab="Predicted Population", xlab="Observed Population",
                   legend = "right", legend.title=element_text("Count",size=22) ,
                   font.legend=c(20),
                   legend.show=F,
                  font.label = list(size = 20, face = "bold", color ="red"),
                  font.x = c(20),
                  font.y = c(20),
                  font.main=c(16),
                  font.xtickslab =c(20),
                  font.ytickslab =c(20),
                  #xticks.by = TRUE,
                  xlim= c(0, max(df1$x)),
                   xtickslab.rt = 45, ytickslab.rt = 45)

plot_cup2

################################################
# District 
#### - BHM
lgdt2 <- cu2 %>% group_by(Dist_Name)%>% drop_na(mean) %>%
  summarise(Obs = sum(POPN, na.rm=T),
            Pred = round(sum(mean, na.rm=T)))
lgdt2 <- lgdt2[,c("Obs", "Pred")]
df2 = data.frame(lgdt2); colnames(df2) = c("x","y")

cup2 <-ggplot(data=df2,aes(x,y)) + 
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
  scale_fill_continuous(type="viridis") +
  geom_smooth(method=lm,linetype=2,colour="red",se=F) + 
  scale_y_continuous(breaks=seq(0,max(df2$y, na.rm=T),50000))+
  guides(alpha="none") +
  geom_point() + commonTheme

round(cor(df2$x,df2$y),2)

plot_dst2<- ggpar(cup2, ylab="Predicted Population", xlab="Observed Population",
                  legend = "right", legend.title=element_text("Density",size=22) ,
                  font.legend=c(20),
                  font.label = list(size = 20, face = "bold", color ="red"),
                  font.x = c(20),
                  font.y = c(20),
                  font.main=c(16),
                  font.xtickslab =c(20),
                  font.ytickslab =c(20),
                  #xticks.by = TRUE,
                  xlim= c(0, max(df2$x)),
                  xtickslab.rt = 45, ytickslab.rt = 45)

plot_dst2

## TSBHM 
lgdt <- cu1 %>% group_by(Dist_Name)%>% drop_na(mean) %>%
  summarise(Obs = sum(POPN, na.rm=T),
            Pred = round(sum(mean, na.rm=T)))
lgdt1 <- lgdt[,c("Obs", "Pred")]
df1 = data.frame(lgdt1); colnames(df1) = c("x","y")
cup1 <-ggplot(data=df1,aes(x,y)) + # swapped 
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
  scale_fill_continuous(type="viridis") +
  geom_smooth(method=lm,linetype=2,colour="red",se=F) + 
  scale_y_continuous(breaks=seq(0, max(df1$y),50000))+
  # scale_x_continuous(breaks=seq(0,max(df1$x, na.rm=T),50000))+
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
                  xlim= c(0, max(df1$x)),
                  xtickslab.rt = 45, ytickslab.rt = 45)

plot_dst1

cor(df1$x, df1$y)
ggarrange(plot_cup1, plot_dst1+ rremove("x.text"), 
          labels = c("(A)", "(B)"),
          nrow=2)



#------------------------------------------------------------------------------

ggdensity(wdata, x = "weight",
          add = "mean", rug = TRUE,
          color = "sex", fill = "sex",
          palette = c("#00AFBB", "#E7B800"))

#---save
#write.csv(cu.data, file=paste0(results_path, "/updated2/CU_data_with_estimates.csv"))

###
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



#factor(sample(x = rep(1:5, each = floor(32090 / 5)),  # Sample IDs for training data
 #             size = 32090))
#
#################################################################################
#####-------------------Cross validation
####################################################################################

cross_val <- function(dat, mod, mesh, spde,
                      A, shp, formula, k_folds)
{
  set.seed(9444506)                                     
  N <- nrow(dat)
  # (cv <- cross_val(datt, mod1c, fit11c, mesh, spde)
  
  dat = datt
  ######
  ind_train <- factor(sample(x = rep(1:k_folds, each = floor(N/ k_folds)),  # Sample IDs for training data
                             size = N, rep=T))
  
  table(as.numeric(ind_train)) 
  dat$k_fold <- as.numeric(ind_train)
  coords <- cbind(dat$lon, dat$lat)
  
  
  k_uniq <-sort(unique(dat$k_fold))
  metrics_cv <- matrix(0, nrow=length(k_uniq), ncol=5)#5 metrics for each fold
  
  cvs <- list()
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
    ###################################
    
    ##==========================
    
    ###---Rerun INLA for model test prediction
    model <-inla(formula, #the formula
                 data=inla.stack.data(stk_train,spde=spde),  #the data stack
                 family= 'gamma',   #which family the data comes from
                 control.predictor=list(A=inla.stack.A(stk_train),compute=TRUE),  #compute gives you the marginals of the linear predictor
                 control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                 verbose = FALSE) #can include verbose=TRUE to see the log of the model runs
    summary(model)
    
    
    
    
    #######################----Extract Spatial Random effects
    sfield_nodes_mean <- mod$summary.random$spatial.field['mean']
    field_mean <- (A%*% as.data.frame(sfield_nodes_mean)[, 1])
    
    sfield_nodesL <- mod$summary.random$spatial.field['0.025quant']
    fieldL <- (A%*% as.data.frame(sfield_nodesL)[, 1])
    
    sfield_nodesU<- mod$summary.random$spatial.field['0.975quant']
    fieldU <- (A%*% as.data.frame(sfield_nodesU)[, 1])
    
    
    
    ###----Extract settlement type random effects
    set <- set_t(dat, mod$summary.random$set_typ$mean)
    setL <- set_t(dat, mod$summary.random$set_typ$`0.025quant`)
    setU <- set_t(dat, mod$summary.random$set_typ$`0.975quant`)
    
    
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
      #rnorm(nrow(test),0, 1/model$summary.hyperpar$mean[1]) + #
      field_mean[-train_ind,1]
    
    dens_ht <- exp(fixed)
    sum(pop_ht <- dens_ht*test$bld_pred)
    
    
    ###
    
    par(mfrow=c(2,2))
    plot(test$POPN, pop_ht, xlab = "Observed", 
         ylab = "Predicted", col=c('dark green','orange'),
         pch=c(16,16), cex.axis=1.5)
    abline(0,1)
    legend("topleft", c("Observed", "Predicted"), col=c("dark green", "orange"), pch=c(16,16),
           bty="n", cex=1.5) 
    
    
    
    (met <- mod_metrics(test$POPN,  
                          pop_ht, pop_htU,  pop_htL, mod))
    
    
    
    
    metrics_cv[i, ] <- as.vector(unlist(met)) 
    
    cvs[[i]] <- data.frame(test=test$POPN, pred=pop_ht, fold=rep(k_uniq[i], length(as.vector(pop_ht))))
  }
  
  
  stat = data.frame(metrics =c("MAE", "RMSE", "BIAS","IMPRECISION","Accuracy"))
  values = data.frame(metrics_cv)
  rownames(values) = paste0("Fold_", 1:k_folds, sep="") 
  
  metrics = bind_cols(stat, values)
  
  return(list(metrics= metrics, data = cvs))
}

#print(nm)
k_folds=5

(cv1 <- cross_val(datt, mod1a, mesh, spde,#--TSBHM
                 A, shp, form1a, k_folds))

datt4 <- datt
datt4$dens2 <- datt$dens3
(cv2 <- cross_val(datt4, mod1b, mesh, spde,#--BHM
                  A, shp, form1b, k_folds))




cv1$metrics
cv2$metrics
write.csv(cv, paste0(out_path, "/cross_validation.csv"))



#mean(c(11781559, 11998136, 11852921, 11745906))

# Basic scatter plot.
dt1 <- data.frame(test = datt$POPN, pred=fit11a, fold=rep(0, nrow(datt)))
dtt1 <- bind_rows(dt1, cv1$data[[1]], cv1$data[[2]],cv1$data[[3]],
                 cv1$data[[4]],cv1$data[[5]])

dtt1$fold1 <- factor(dtt1$fold)
levels(dtt1$fold1) <- paste0("fold", 0:5, sep="")
table(dtt1$fold1)


##----Save workspace
save.image(paste0(results_path, "/workspace_gg_model.Rdata"))
