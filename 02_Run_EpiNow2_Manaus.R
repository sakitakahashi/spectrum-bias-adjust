library(dplyr)
library(tidyr)
library(tidyverse)
library(stringr)
library(EpiNow2)
library(lubridate)
library(extraDistr)
library(readxl)

## Read in the data
read_excel("Data/TimeSeries_data_Manaus_hospitalizations.xlsx", sheet=1) %>%
	mutate(date=as.Date(date)) %>%
	arrange(date) %>%
	filter(date >= ymd("2020-03-01")) %>%
	complete(date=seq.Date(min(date), max(date), by="day"), fill=list(hosp=0)) -> dat_ts_Manaus

## Set up example generation time
generation_time <- get_generation_time(disease="SARS-CoV-2", source="ganyani")

###############################################################
## [1A] Set delays between symptom onset and HOSPITAL report ##
## Parameters from Bi et al (2020):							             ##
## https://doi.org/10.1016/S1473-3099(20)30287-5			       ##
###############################################################

## Parametrize
Bi_symptomOnset_to_hosp_param_1 <- 1.23
Bi_symptomOnset_to_hosp_param_2 <- 0.79

lag_symptomOnset_to_hosp_Bi <- bootstrapped_dist_fit(rgamma(1000, Bi_symptomOnset_to_hosp_param_1, Bi_symptomOnset_to_hosp_param_2), dist="lognormal", max_value=14)

###############################################################
## [1B] Set delays between symptom onset and HOSPITAL report ##
## Parameters from Linton et al (2020):						           ##
## https://doi.org/10.3390/jcm9020538						             ##
###############################################################

## Parametrize
Linton_symptomOnset_to_hosp_mean <- 9.7
Linton_symptomOnset_to_hosp_SD <- 35.2

## Get lognormal parameters from mean and sd
## Code adapted from: https://github.com/jriou/covid_adjusted_cfr/blob/master/run_models/run_model16_spain.R
get_par_lnorm <- function(m,s) {
	mu <- log(m) - 1/2*log((s/m)^2+1)
	sigma <- sqrt(log((s/m)^2+1))
	return(list(mu=mu, sigma=sigma))
}

## Get lognormal parameters from mean and sd
Linton_symptomOnset_to_hosp_param_1 <- get_par_lnorm(Linton_symptomOnset_to_hosp_mean, Linton_symptomOnset_to_hosp_SD)$mu
Linton_symptomOnset_to_hosp_param_2 <- get_par_lnorm(Linton_symptomOnset_to_hosp_mean, Linton_symptomOnset_to_hosp_SD)$sigma

lag_symptomOnset_to_hosp_Linton <- bootstrapped_dist_fit(rlnorm(1000, Linton_symptomOnset_to_hosp_param_1, Linton_symptomOnset_to_hosp_param_2), dist="lognormal", max_value=70)

################################################################################

############################################################
## Run model without Rt estimation (just backcalculation) ##
############################################################

## Clean up the data
dat_ts_Manaus %>%
	rename(confirm = hosp) %>%
	select(date, confirm) -> HOSP_DATA

## Fit
backcalc_using_reported_hosp_Bi <- estimate_infections(
	reported_cases=HOSP_DATA,
	generation_time=generation_time,
	family="negbin",
	CrIs=c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975),
	week_effect=TRUE,
	rt_prior=list(),
	delays=list(lag_symptomOnset_to_hosp_Bi),
	return_fit=TRUE,
	verbose=TRUE,
	use_breakpoints=FALSE,
	stan_args=list(warmup=2000, chains=4, cores=1))

backcalc_using_reported_hosp_Linton <- estimate_infections(
	reported_cases=HOSP_DATA,
	generation_time=generation_time,
	family="negbin",
	CrIs=c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975),
	week_effect=TRUE,
	rt_prior=list(),
	delays=list(lag_symptomOnset_to_hosp_Linton),
	return_fit=TRUE,
	verbose=TRUE,
	use_breakpoints=FALSE,
	stan_args=list(warmup=2000, chains=4, cores=1))

## Save the plot objects
plot_hosp_Bi <- plot_estimates(estimate=backcalc_using_reported_hosp_Bi$summarised[variable=="infections"], reported=HOSP_DATA)

plot_hosp_Linton <- plot_estimates(estimate=backcalc_using_reported_hosp_Linton$summarised[variable=="infections"], reported=HOSP_DATA)

## Save
save(backcalc_using_reported_hosp_Bi, backcalc_using_reported_hosp_Linton, plot_hosp_Bi, plot_hosp_Linton, file="Results/EpiNow2_fit_Manaus.RData")
