library(dplyr)
library(tidyr)
library(tidyverse)
library(stringr)
library(EpiNow2)
library(lubridate)
library(extraDistr)
library(readxl)

## Read in the data
load("Data/TimeSeries_data_Japan.RData")

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

## Clean up the data
dat_ts_Japan %>%
	filter(date <= ymd("2021-01-15")) %>%
	droplevels() %>%
	rename(confirm = pcr_positive) %>%
	select(date, confirm) -> CASE_DATA

dat_ts_Japan %>%
	filter(date <= ymd("2021-01-15")) %>%
	droplevels() %>%
	rename(confirm = deaths_delta) %>%
	select(date, confirm) -> DEATH_DATA

## Fit
backcalc_using_reported_cases <- estimate_infections(
	reported_cases=CASE_DATA,
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
	reported_cases=DEATH_DATA,
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
plot_cases <- plot_estimates(estimate=backcalc_using_reported_cases$summarised[variable=="infections"], reported=CASE_DATA, ylab="Cases")

plot_deaths <- plot_estimates(estimate=backcalc_using_reported_deaths$summarised[variable=="infections"], reported=DEATH_DATA, ylab="Deaths")

## Save
save.image(file="Results/EpiNow2_fit_Japan.RData")
