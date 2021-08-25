library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(readxl)
library(rstan)
library(loo)
library(gridExtra)
library(scales)
library(doParallel)

## How many processors?
CORES <- 1
registerDoParallel(CORES)

## Set input parameters
input_params <- data.frame(FILENAME_TO_READ=rep(NA,2), FILENAME_TO_WRITE=NA, ALGO=NA)

## Filenames to read & to write
input_params$FILENAME_TO_READ[1] <- "Data/Kinetics_data_bivariate_2groups.RData"
input_params$FILENAME_TO_READ[2] <- "Data/Kinetics_data_bivariate_3groups.RData"

input_params$FILENAME_TO_WRITE[1] <- "Results/Stan_fit_bivariate_2groups.RData"
input_params$FILENAME_TO_WRITE[2] <- "Results/Stan_fit_bivariate_3groups.RData"

## How many disease severity groups to model?
input_params$GROUPS[1] <- 2
input_params$GROUPS[2] <- 3

## Which Stan algorithm to use?
input_params$ALGO[1] <- "HMC"
input_params$ALGO[2] <- "HMC"

## Loop over all 
for(i in 1:nrow(input_params)) {
  
  ## Load the data & compiled code
  load(file=input_params$FILENAME_TO_READ[i])
  
  ## Run Stan iterations
  if(input_params$GROUPS[i]==2) {
      
    fit_Stan <- stan(
      file="Code/bivariate_mixture_log_model_2groups.stan",
      data=data_Stan,
      algorithm=input_params$ALGO[i], seed="12345", cores=CORES, chains=CORES, iter=50000, thin=50, init=0)
    
  }
  
  if(GROUPS[i]==3) {
    
    fit_Stan <- stan(file="Code/bivariate_mixture_log_model_3groups.stan",
      data=data_Stan,
      algorithm=input_params$ALGO[i], seed="12345", cores=CORES, chains=CORES, iter=50000, thin=50, init=0)
  
  }

  ## Save
  save(fit_Stan, file=input_params$FILENAME_TO_WRITE[i])
  
  rm(fit_Stan)

}
