library(dplyr)
library(tidyr)
library(tidyverse)
library(stringr)
library(EpiNow2)
library(lubridate)
library(extraDistr)

## Read in the data
load("Data/TimeSeries_data_US_NYT.RData")

## Get the state names
STATES <- levels(data_NYT_cleaned$state)

## Set up example generation time
generation_time <- get_generation_time(disease="SARS-CoV-2", source="ganyani")

#########################################################
## [1] Set delay between symptom onset and CASE report ##
#########################################################

## Parametrize
Bi_symptomOnset_to_caseReport_param_1 <- 2.12
Bi_symptomOnset_to_caseReport_param_2 <- 0.39

lag_symptomOnset_to_caseReport_Bi <- bootstrapped_dist_fit(rgamma(1000, Bi_symptomOnset_to_caseReport_param_1, Bi_symptomOnset_to_caseReport_param_2), dist="lognormal", max_value=14)

###########################################################
## [2] Set delays between symptom onset and DEATH report ##
###########################################################

## Parametrize
Linton_symptomOnset_to_death_mean <- 20.2
Linton_symptomOnset_to_death_SD <- 11.6

## Get lognormal parameters from mean and sd
## Code adapted from: https://github.com/jriou/covid_adjusted_cfr/blob/master/run_models/run_model16_spain.R
get_par_lnorm <- function(m,s) {
	mu <- log(m) - 1/2*log((s/m)^2+1)
	sigma <- sqrt(log((s/m)^2+1))
	return(list(mu=mu, sigma=sigma))
}

## Get lognormal parameters from mean and sd
Linton_symptomOnset_to_death_param_1 <- get_par_lnorm(Linton_symptomOnset_to_death_mean, Linton_symptomOnset_to_death_SD)$mu
Linton_symptomOnset_to_death_param_2 <- get_par_lnorm(Linton_symptomOnset_to_death_mean, Linton_symptomOnset_to_death_SD)$sigma

lag_symptomOnset_to_death_Linton <- bootstrapped_dist_fit(rlnorm(1000, Linton_symptomOnset_to_death_param_1, Linton_symptomOnset_to_death_param_2), dist="lognormal", max_value=70)

################################################################################

############################################################
## Run model without Rt estimation (just backcalculation) ##
############################################################

for(i in 1:52) {
	
	## Clean up the data
	data_NYT_cleaned %>%
		filter(state==STATES[i]) %>%
		# filter(date <= ymd("2020-10-14")) %>%
		droplevels() %>%
		rename(confirm = new_cases) %>%
		select(date, confirm) -> CASE_DATA_FOR_STATE
	
	data_NYT_cleaned %>%
		filter(state==STATES[i]) %>%
		# filter(date <= ymd("2020-10-14")) %>%
		droplevels() %>%
		rename(confirm = new_deaths) %>%
		select(date, confirm) -> DEATH_DATA_FOR_STATE
	
	## Fit
	backcalc_using_reported_cases <- estimate_infections(
		reported_cases=CASE_DATA_FOR_STATE,
		generation_time=generation_time,
		family="negbin",
		CrIs=c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975),
		week_effect=TRUE,
		rt_prior=list(),
		delays=list(lag_symptomOnset_to_caseReport_Bi),
		return_fit=TRUE,
		verbose=TRUE,
		use_breakpoints=FALSE,
		stan_args=list(warmup=2000, chains=4, cores=1))
	
	backcalc_using_reported_deaths <- estimate_infections(
		reported_cases=DEATH_DATA_FOR_STATE,
		generation_time=generation_time,
		family="negbin",
		CrIs=c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975),
		week_effect=TRUE,
		rt_prior=list(),
		delays=list(lag_symptomOnset_to_death_Linton),
		return_fit=TRUE,
		verbose=TRUE,
		use_breakpoints=FALSE,
		stan_args=list(warmup=2000, chains=4, cores=1))
	
	## Save the plot objects
	plot_cases <- plot_estimates(estimate=backcalc_using_reported_cases$summarised[variable=="infections"], reported=CASE_DATA_FOR_STATE, ylab="Cases")
	
	plot_deaths <- plot_estimates(estimate=backcalc_using_reported_deaths$summarised[variable=="infections"], reported=DEATH_DATA_FOR_STATE, ylab="Deaths")
	
	## Save
	save.image(file=paste0("Results/EpiNow2_fit_US_", STATES[i], ".RData"))
	
	## Clean up
	print(i)
	rm(CASE_DATA_FOR_STATE)
	rm(DEATH_DATA_FOR_STATE)
	
}

################################################################################

## Load in a state EpiNow2 run - cases & deaths
load(paste0("Results/EpiNow2_fit_US_", STATES[1], ".RData"))

backcalc_using_reported_cases["summarised"] %>%
  as.data.frame() %>%
  filter(summarised.variable=="infections") %>%
  droplevels() %>%
  mutate(summarised.date = ymd(summarised.date)) %>%
  mutate(state=STATES[1]) -> STORAGE_CASES_EpiNow2

CASE_DATA_FOR_STATE %>%
  as.data.frame() %>%
  mutate(state=STATES[1]) -> REPORTED_CASES

backcalc_using_reported_deaths["summarised"] %>%
  as.data.frame() %>%
  filter(summarised.variable=="infections") %>%
  droplevels() %>%
  mutate(summarised.date = ymd(summarised.date)) %>%
  mutate(state=STATES[1]) -> STORAGE_DEATHS_EpiNow2

DEATH_DATA_FOR_STATE %>%
  as.data.frame() %>%
  mutate(state=STATES[1]) -> REPORTED_DEATHS

## Load in all the state EpiNow2 runs
for(i in 2:length(STATES)) {
  
  load(paste0("Results/EpiNow2_fit_US_", STATES[i], ".RData"))
  
  backcalc_using_reported_cases["summarised"] %>%
    as.data.frame() %>%
    filter(summarised.variable=="infections") %>%
    droplevels() %>%
    mutate(summarised.date = ymd(summarised.date)) %>%
    mutate(state=STATES[i]) -> tmp_cases
  
  STORAGE_CASES_EpiNow2 <- bind_rows(STORAGE_CASES_EpiNow2, tmp_cases)
  
  CASE_DATA_FOR_STATE %>%
    as.data.frame() %>%
    mutate(state=STATES[i]) -> rep_cases
  
  REPORTED_CASES <- bind_rows(REPORTED_CASES, rep_cases)
  
  backcalc_using_reported_deaths["summarised"] %>%
    as.data.frame() %>%
    filter(summarised.variable=="infections") %>%
    droplevels() %>%
    mutate(summarised.date = ymd(summarised.date)) %>%
    mutate(state=STATES[i]) -> tmp_deaths
  
  STORAGE_DEATHS_EpiNow2 <- bind_rows(STORAGE_DEATHS_EpiNow2, tmp_deaths)
  
  DEATH_DATA_FOR_STATE %>%
    as.data.frame() %>%
    mutate(state=STATES[i]) -> rep_deaths
  
  REPORTED_DEATHS <- bind_rows(REPORTED_DEATHS, rep_deaths)
  
  print(i)
  
}

## Save everything
save(STORAGE_CASES_EpiNow2, STORAGE_DEATHS_EpiNow2, REPORTED_CASES, REPORTED_DEATHS, file="Results/EpiNow2_fit_US.RData")
