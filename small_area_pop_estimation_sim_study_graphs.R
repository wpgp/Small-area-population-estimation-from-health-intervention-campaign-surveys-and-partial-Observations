
#----------------------------------------------------------------------------------------------------
# Posterior Exploration of the simulated data
#--------------------------------------------------------------------------------------------------------
# Load Results and explore the posterior samples/data
#----------------------------------------------------------------------------------------------------------
# load required packages (Note that most of these packages are already loaded during simulation)
library(INLA); library(raster); library(maptools)
library(gtools); library(sp); library(spdep); library(rgdal)
library(fields); library(mvtnorm);  library(geoR)
library(actuar);library(viridisLite);require(grid);require(gridExtra)
require(lattice);require(tidyverse);require(MASS);library(tmap)
library(tmaptools);library(sf)
library(cartography) # mapping dedicated package
library(OpenStreetMap)


# Set the directory where the simulated data are stored so they are easily retrieved
path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP000008_UNFPA_PNG/Working/Chris/reviews/sim_study2" # please change to your own working directory
out_path <- paste0(path, "/outputs")
file_path <- paste0(path, "/file")
setwd(out_path)


# samples
pop.cover <- c(100, 80, 60, 40, 20) # survey coverage 
bld.cover <- c(100, 95, 90, 85, 80, 75, 70, 65)# proportion of settlement data observed

metric <- c("mae", "rmse", "bias", "corr") # model fit metrics to calculate
method <- c("onestep", "twostep") # onestep - BHM; twostep = TSBHM

n.pop.cover <- length(pop.cover)
n.bld.cover <- length(bld.cover)
n.metric <- length(metric)
n.method <- length(method)

##---build the dataframe for the metrics
n.size <- n.pop.cover*n.bld.cover*n.metric*n.method
dim(dat.met <- data.frame(expand.grid(method=method,
                                      bld_cover=bld.cover,
                                      pop_cover=pop.cover,
                                      metric=metric)))

###
# Extract the model fit metrics
dat_met1 <- list()
dat_met2 <- list()
for(j in 1:n.pop.cover)
{
  pathp <- paste0(out_path,"/outputs_for_", pop.cover[j],"%","_pop_count")
  for(k in 1:n.bld.cover)
  {
    pathb <- paste0(pathp,"/", bld.cover[k], "%","_bldg_count")
    met0 <- read.csv(paste0(pathb, "/fit_metrics_0.csv"))  # bhm 
    met1 <- read.csv(paste0(pathb, "/fit_metrics_1.csv")) # tsbhm
    met0[1] <- bld.cover[k]
    met1[1] <- bld.cover[k]
    met0 <- c(pop.cover[j], met0)
    met1 <- c(pop.cover[j], met1)
    dat_met1[[k]] = rbind(met0, met1) 
  }
  dat_met2[[j]] = dat_met1
} 
#dat_met2
unnest_list <- unlist(dat_met2, recursive = FALSE)  #--unnest the list
str(unnest_list)
dim(metrics <- as.data.frame(matrix(unlist(do.call(rbind, unnest_list)),
                                    nrow=80, ncol=6)))
##
names(metrics) <- c("pop_cover", "bld_cover", "mae","rmse", "bias", "corr") #-rename columns
metrics$method <- rep(c("onestep", "twosteps"),nrow(metrics)/2)#--add 'method' col

head(metrics)
write.csv(metrics, "combined_fit_metrics.csv", row.names=FALSE)

# Convert to long format for plotting 
require(reshape2)
require(ggpubr)
dim(met_long <- melt(metrics, id.vars=c("pop_cover","bld_cover", "method"),
                     value.name="estimate", variable.name = "metric"))

met_long$method = factor(met_long$method)
met_long$pop_cover = factor(met_long$pop_cover)
head(met_long)
table(met_long$metric)
#write.csv(met_long, "combined_fit_metrics_long.csv", row.names=FALSE)

##---Variable Recode
variable_names <- list(
  "mae" = "Mean  Absolute \n Error (MAE)" ,
  "rmse" = "Root Mean Square \n Error (RMSE)",
  "bias" = "Absolute Bias",
  "corr" = "Correlation \n Coefficient"
)

levels(met_long$method) <- c("BHM","TSBHM") # rename 
grp <- levels(met_long$method)

variable_labeller2 <- function(variable,value){
  if (variable=='metric') {
    return(variable_names[value])
  } else {
    return(grp)
  }
}



# Group plot
####
var_names <- c(
  "100",
  "80",
  "60",
  "40" ,
  "20"
)

met_long$metric2 <- factor(met_long$metric, levels=var_names)

table(met_long$metric)

variable_names <- list(
  "mae" = "Mean Absolute \n Error",
  "rmse" = "Root Mean Square \n Error",
  "bias" = "Absolute Bias",
  "corr" = "Correlation \n Coefficient"
)

grp <- levels(met_long$method)

variable_labeller <- function(variable,value){
  if (variable=='metric2') {
    return(variable_names[value])
  } else {
    return(grp)
  }
}
###


####============------------------------
#---------------Make Scatter plots for pop counts------------

Var2Include <- c("lon", "lat", "prov2_ID","bld", "pop",
                 "dens", "popm", "bldm","mean_dens_hat", 
                 "mean", "lower", "upper",
                 "pop_cover", "bld_cover", "method")

dat_cua<- list()
dat_cub <- list()
for(j in 1:n.pop.cover)
{
  # set initial directory 
  pathp <- paste0(out_path,"/outputs_for_", pop.cover[j],"%","_pop_count")
  for(k in 1:n.bld.cover)
  {
    
    # set the nested directory
    pathb <- paste0(pathp,"/", bld.cover[k], "%","_bldg_count")
    
    
    # load cu level data
    cu0 <- read.csv(paste0(pathb, "/CU_estimates_0.csv")) # bhm
    cu1 <- read.csv(paste0(pathb, "/CU_estimates_1.csv")) #tsbhm
    
    # Add method variable 
    cu0$method <- rep("BHM", nrow(cu0))
    cu1$method <- rep("TSBHM", nrow(cu1))
    
    # 
    cu0$bld_cover <- rep(bld.cover[k], nrow(cu0))
    cu0$pop_cover <- rep(pop.cover[j], nrow(cu0))
    
    cu1$bld_cover <- rep(bld.cover[k], nrow(cu1))
    cu1$pop_cover <- rep(pop.cover[j], nrow(cu1))
    
    cu0 <- cu0[,Var2Include]
    cu1 <- cu1[,Var2Include]
    
    dat_cua[[k]] = rbind(cu0,cu1)
  }
  dat_cub[[j]] = dat_cua
} 

# unlist and unest the nested lists
unnest_cu <- unlist(dat_cub, recursive = FALSE)  #--unnest the list
str(unnest_cu)
dim(cu_dat <- as.data.frame(matrix(unlist(do.call(rbind, unnest_cu)),
                                   nrow=64200*8*5, ncol=15))) # 2,568,00 by 14 

colnames(cu_dat) <- Var2Include # variable name

#write.csv(cu_dat, "cu_posterior_samples.csv")

require(ggpubr)


#####################
###
var_names <- c(
  "100",
  "80",
  "60",
  "40" ,
  "20"
)

cu_dat$pop_cover2 <- factor(cu_dat$pop_cover, levels=var_names)


variable_names <- list(
  "100" = "100% \n Survey Coverage",
  "80" = "80% \n Survey Coverage",
  "60" = "60% \n Survey Coverage",
  "40" = "40% \n Survey Coverage",
  "20" = "20% \n Survey Coverage" 
)

grp <- levels(cu_dat$method)

variable_labeller <- function(variable,value){
  if (variable=='pop_cover2') {
    return(variable_names[value])
  } else {
    return(grp)
  }
}


#--------------------------plot all survey coverage and all satellite coverage props together
# FIGURE S6 (Supplemental)

fdat100 <- CU_dat %>% filter(bld_cover==100)
class(fdat100$method)
names(fdat100)
p100 <- fdat100 %>%
  ggplot(aes(x=pop, y=mean))+
  geom_point()+
  geom_errorbar(aes(ymin = lower, ymax = upper,colour=pop_cover), size=1,width = 2)+
  geom_smooth(aes(colour=pop_cover),method="lm", se=T)+
  theme_bw()+
  theme(strip.text = element_text(size = 15),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  facet_wrap(~method)
p100

plot100 <- ggpar(p100, xlab="Observed Population(Count)", ylab="Predicted Population(Count)",
                 legend = "right", legend.title = "Survey \n Coverage (%)",
                 font.legend=c(20),
                 palette = c("jco"),
                 yscale = c("none"),
                 font.label = list(size = 18, face = "bold", color ="red"),
                 font.x = c(22),
                 font.y = c(20),
                 font.main=c(14),
                 font.xtickslab =c(18),
                 font.ytickslab =c(20),
                 xtickslab.rt = 45, ytickslab.rt = 45)

plot100


##-------------------------------------------------------

##---------------------------------------------------------------------------------
##--recode survey coverage levels
levels(fdat$pop_cover2) <- c("Survey \n Coverage:100%", 
                             "Survey \n Coverage:80%",
                             "Survey \n Coverage:60%", 
                             "Survey \n Coverage:40%", 
                             "Survey \n Coverage:20%")

##--Recode Satellite Observations levels
var_namesb <- c(
  "100",
  "95",
  "90",
  "85" ,
  "80",
  "75",
  "70" ,
  "65"
)
fdat$bld_cover2 <- factor(fdat$bld_cover, levels=var_namesb)
table(fdat$bld_cover2)

levels(fdat$bld_cover2)  <- c(
  "Satellite \n Coverage:100%",
  "Satellite \n Coverage:95%",
  "Satellite \n Coverage:90%",
  "Satellite \n Coverage:85%",
  "Satellite \n Coverage:80%", 
  "Satellite \n Coverage:75%",
  "Satellite \n Coverage:70%",
  "Satellite \n Coverage:65%"
)

names(cu_dat)
nat_dat1 <- cu_dat %>% dplyr::select(mean, lower, upper, pop_cover, bld_cover,
                                     method, pop_cover2) %>%
  summarise()
#----------------------------------------------------------------
#  Extract the notal national estimates obtained from the various data situations 
# This comes with uncertainty estimates which could not be obtained by simply agggregating
# the cu data
#---------------National estimates------------

Var2Include <- c("lon", "lat", "prov2_ID","bld", "pop",
                 "dens", "popm", "bldm","mean_dens_hat", 
                 "mean", "lower", "upper",
                 "pop_cover", "bld_cover", "method")
dat_nata<- list()
dat_natb <- list()
for(j in 1:n.pop.cover)
{
  
  # specify outer data path
  pathp <- paste0(out_path,"/outputs_for_", pop.cover[j],"%","_pop_count")
  for(k in 1:n.bld.cover)
  {
    
    # specify inner data path
    pathb <- paste0(pathp,"/", bld.cover[k], "%","_bldg_count")
    
    # load the national estimates popsterior samples 
    nat0 <- read.csv(paste0(pathb, "/national_estimates_0.csv"))  #bhm
    nat1 <- read.csv(paste0(pathb, "/national_estimates_1.csv")) # tsbhm
    
    # Add method variable 
    nat0$method <- rep("BHM", nrow(nat0))
    nat1$method <- rep("TSBHM", nrow(nat1))
    
    # 
    nat0$bld_cover <- rep(bld.cover[k], nrow(nat0))
    nat0$pop_cover <- rep(pop.cover[j], nrow(nat0))
    
    nat1$bld_cover <- rep(bld.cover[k], nrow(nat1))
    nat1$pop_cover <- rep(pop.cover[j], nrow(nat1))
    
    
    dat_nata[[k]] = rbind(nat0,nat1)
  }
  dat_natb[[j]] = dat_nata
} 

str(dat_natb)
# unlist and unest the nested lists
unnest_nat <- unlist(dat_natb, recursive = FALSE)  #--unnest the list
str(unnest_nat)
dim(nat_dat <- as.data.frame(matrix(unlist(do.call(rbind, unnest_nat)),
                                    nrow=8*8*5, ncol=5))) # 2,568,00 by 14 

# variables recode
colnames(nat_dat) <- c("measure", "estimate", "method", "bld_cover", "pop_cover")

class(nat_dat$measure <- factor(nat_dat$measure))
levels(nat_dat$measure) <- c("lower", "median", "total", "upper")


# subset  and combine key datasets
dim(nat_mean <- nat_dat[nat_dat$measure=="total",])
dim(nat_lower <- nat_dat[nat_dat$measure=="lower",])
dim(nat_upper <- nat_dat[nat_dat$measure=="upper",])
nat_mean$lower <- round(as.numeric(nat_lower$estimate))
nat_mean$upper <- round(as.numeric(nat_upper$estimate))
nat_mean$estimate <- round(as.numeric((nat_mean$estimate)))



###---plots-----------------------
var_names <- c(
  "100",
  "80",
  "60",
  "40" ,
  "20"
)

nat_mean$pop_cover2 <- factor(nat_mean$pop_cover, levels=var_names)


variable_names <- list(
  "100" = "100% \n Survey Coverage",
  "80" = "80% \n Survey Coverage",
  "60" = "60% \n Survey Coverage",
  "40" = "40% \n Survey Coverage",
  "20" = "20% \n Survey Coverage" 
)

ngrp <- levels(nat_mean$method)

variable_labeller <- function(variable,value){
  if (variable=='pop_cover2') {
    return(variable_names[value])
  } else {
    return(ngrp)
  }
}


levels(nat_mean$pop_cover2) <- c(
  "Survey Coverage:\n 100%",
  "Survey Coverage: \n 80%",
  "Survey Coverage: \n 60%",
  "Survey Coverage: \n 40%",
  "Survey Coverage: \n 20%" 
)

# Create a simple example dataset

#bld_cover



# -----------------------------------------------------
nat_mean100 <- nfdat_mean[nat_mean$pop_cover== "100",]
nat_mean100$mean2 <- nat_mean100$mean/1000000


# -----------------------------------------------------

nat_mean60 <- nat_mean[nat_mean$pop_cover== "60",]
nat_mean60$mean2 <- nat_mean60$estimate/1000000
nat_mean60$upper2 <- nat_mean60$upper/1000000
nat_mean60$lower2 <- nat_mean60$lower/1000000


####------------------------------------------------------
# N                          figure 3
#---------------------------------------------------------------
##----- RMSE # Fig 3A
dim(rmse <- met_long[met_long$metric=="rmse",])

plot_rmse <- ggline(rmse, x = "bld_cover", y = "estimate",
                    error.plot = "estimate",
                    facet.by = "method",
                    panel.labs= list(method=c("BHM", "TSBHM")),
                    panel.labs.font.x = list(size=20),
                    color = "pop_cover",
                    point.size=1.5,
                    #linetype = "pop_cover",
                    size=1.4)
rrmse <-  ggpar(plot_rmse, xlab="Settlement proportion observed (%)", ylab="Root mean square error (RMSE)",
                legend = "top", legend.title = "Survey \n Coverage (%)",size=22,
                font.legend=c(18),
                # palette = c("#00AFBB", "#E7B800", "#FC4E07", "#0D0887FF", "#993333"),
                palette = "jco",
                #colour = "bld_cover",
                #shape= "bld_cover",
                font.label = list(size = 15, face = "bold", color ="red"),
                font.x = c(22),
                font.y = c(20),
                font.main=c(20),
                font.xtickslab =c(16),
                font.ytickslab =c(20),
                # orientation = "reverse",
                xtickslab.rt = 45, ytickslab.rt = 45)
rrmse 

##----- CORR #fig3B
dim(corr <- met_long[met_long$metric=="corr",])
#plot_corr <- ggplot(corr, aes(x=bld_cover, y=estimate, color = pop_cover))+
#geom_point()+
#geom_line()+
#facet_wrap(~method)
plot_corr <- ggline(corr, x = "bld_cover", y = "estimate",
                    error.plot = "estimate",
                    facet.by = "method",
                    panel.labs= list(method=c("BHM", "TSBHM")),
                    panel.labs.font.x = list(size=20),
                    color = "pop_cover",
                    point.size=1.5,
                    #linetype = "pop_cover",
                    size=1.4)
rcorr <- ggpar(plot_corr, xlab="Settlement proportion observed (%)", ylab="Correlation coefficient (CC)",
               legend = "top", legend.title = "Survey \n Coverage (%)",size=22,
               font.legend=c(18),
               # palette = c("#00AFBB", "#E7B800", "#FC4E07", "#0D0887FF", "#993333"),
               palette = "jco",
               #colour = "bld_cover",
               #shape= "bld_cover",
               font.label = list(size = 15, face = "bold", color ="red"),
               font.x = c(22),
               font.y = c(20),
               font.main=c(20),
               font.xtickslab =c(16),
               font.ytickslab =c(20),
               # orientation = "reverse",
               xtickslab.rt = 45, ytickslab.rt = 45)
rcorr


# Violin plots  --- # fig 3C
nat_mean100 <- nat_mean %>% filter(pOp_cover==100)
nat_mean100$mean1 <- nat_mean100$estimate/mill
pbx <- ggplot(nat_mean100, aes(x=method, y=mean1, fill=method))+ 
  geom_violin(trim=FALSE)+
  geom_dotplot(binaxis='y', stackdir='center',
               position=position_dodge(1))+
  geom_boxplot(width=0.1, fill="white") +
  geom_hline(yintercept=11625153/1000000 , linetype="dashed", color = "black")+
  theme_bw()

# Scatter plots (sp) with regression linem
rst <- ggpar(pbx,xlab="Method", ylab="Predicted population (millions)",
             legend = "top", legend.title=element_text("Method:"),
             font.legend=c(18),
             palette = c("jco"),
             font.label = list(size = 22, face = "bold", color ="red"),
             font.x = c(22),
             font.y = c(20),
             font.main=c(20),
             font.xtickslab =c(18),
             font.ytickslab =c(20),
             #ylim= c(0, max(nfdat_mean$mean)),
             xtickslab.rt = 45, ytickslab.rt = 45)
rst 


# Multipanel 
ggarrange(rrmse, rcorr, plotp100, rst+ rremove("x.text"), 
          #labels = c("(A)", "(B)", "(C)", "(D)"),
          ncol = 2, nrow = 2)


#-----------------------------------------------------------------------------------
## Calculate Error Rates
#------------------------------------------------------------------------------
dtta <- nat_mean
true <- 11643074


### Subset data for relative error rate calculation 
# BHM
dat_bhm <- nat_mean %>% filter(method == "BHM")%>% 
  filter(pop_cover!=100 & bld_cover!=100)%>% 
  mutate(error = abs(estimate - true)/estimate) # Absolute Error rate


# TSBHM
dat_tsbhm <- nat_mean %>% filter(method == "TSBHM")%>% 
  filter(pop_cover!=100 & bld_cover!=100)%>% 
  mutate(error = abs(estimate - true)/estimate)  # Absolute Error rate


# Relative Error Rate
dat_tsbhm$rer <- dat_tsbhm$error/dat_bhm$error # relative error rate
(dat_tsbhm$prrer <-(1-dat_tsbhm$rer)*100) # percentage reduction in relative error rate

c(min(dat_tsbhm$prrer), max(dat_tsbhm$prrer)) # Range of reduction in relative error rates


# Lollipop plots of reduction in relative error rates caused by the thbhm approach
# Figure 4. 
names(dat_tsbhm)
plot_lp <- dat_tsbhm %>% ggplot(aes(bld_cover, prrer, fill=bld_cover))+ 
  geom_segment(aes(x = bld_cover, xend = bld_cover, y = 0, yend = prrer),
               lwd = 2) +
  geom_point(size =17, pch = 21, bg = "steel blue", col = "black") +
  geom_text(aes(y = round(prrer), label = paste(format(round(prrer)))), size = 6,
            ## make labels left-aligned
            vjust = 1, nudge_y =0.1, nudge_x =0.25,col="black") +
  guides(colour = "none")+ # remove legend
  coord_flip() +
  theme_minimal()+
  theme(strip.text = element_text(size = 18),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14))+
  facet_wrap(~pop_cover2, scales="free", nrow=3)

plp <- ggpar(plot_lp ,xlab="Settlement proportion observed(%)", ylab="Reduction in relative error rate(%)",
             legend = "right", legend.title=element_text("Survey \n coverage(%)"),
             font.legend=c(20),
             palette = c("lancet"),
             font.label = list(size = 22, face = "bold", color ="red"),
             font.x = c(22),
             font.y = c(20),
             font.main=c(20),
             font.xtickslab =c(18),
             font.ytickslab =c(20),
             xtickslab.rt = 45, ytickslab.rt = 45)
plp

