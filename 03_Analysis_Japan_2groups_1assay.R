library(readstata13)
library(dplyr)
library(data.table)
library(ggplot2)
library(tibble)
library(knitr)
library(tidyr)
library(cowplot)
library(stringr)
library(randomForest)
library(RColorBrewer)
library(patchwork)
library(ggrepel)
library(rstan)
library(bayesplot)
library(readr)
library(lubridate)
library(brms)
library(incidental)
library(readxl)
library(GADMTools)
library(scales)
library(jsonlite)
library(tidyverse)
library(stringi)

source("Code/functions_for_gqs.R")

########################
## Load kinetics data ##
########################

## Load each assay
FILENAME_Abbott <- "Results/Stan_fit_Abbott_2groups.RData"
load(FILENAME_Abbott)
fit_Stan_Abbott <- fit_Stan

FILENAME_Roche <- "Results/Stan_fit_Roche_2groups.RData"
load(FILENAME_Roche)
fit_Stan_Roche <- fit_Stan

## Fix parameters
prop_hosp <- 0.05683933

## Get the estimates of specificity
Sp_Abbott <- 0.9963
Sp_Roche <- 0.9980

##############################################################################################
## Calculate severity-weighted & epi curve-weighted sensitivity for each assay, and adjust: ##
## 																							                                            ##
## Japan: using EpiNow2 & deaths															                              ##
##############################################################################################

## Read in merged EpiNow2 data
load("Results/EpiNow2_fit_Japan.RData")

backcalc_using_reported_deaths["summarised"] %>%
	as.data.frame() %>%
	filter(summarised.variable=="infections") %>%
	droplevels() %>%
	mutate(summarised.date = ymd(summarised.date)) %>%
	rename(date = summarised.date) %>%
	rename(onsets = summarised.median) %>%
	select(date, onsets) %>%
	droplevels() -> dat_ts_Japan

ex <- interval(ymd("2020-12-14"), ymd("2020-12-25"))
int_start(ex) + (int_end(ex) - int_start(ex))/2

## Remove 3 weeks
SURVEY_DATE_Japan <- ymd("2020-12-19")-weeks(3)

dat_ts_Japan %>%
	filter(date <= SURVEY_DATE_Japan) %>%
	mutate(days_before_survey = as.double(SURVEY_DATE_Japan - date)) %>%
	mutate(prop_cases = onsets / sum(onsets)) -> dat_ts_Japan_final

## Bring in serosurvey results
dat_sero_Japan_univariate <- read_excel("Data/Seroprevalence_data_Japan.xlsx", sheet=1)

## Get raw prevalence
dat_sero_Japan_univariate$Abbott_raw <- dat_sero_Japan_univariate$Positive_Abbott / dat_sero_Japan_univariate$N_tested
dat_sero_Japan_univariate$Roche_raw <- dat_sero_Japan_univariate$Positive_Roche / dat_sero_Japan_univariate$N_tested

dat_sero_Japan_univariate$Abbott_median <- NA
dat_sero_Japan_univariate$Abbott_lb <- NA
dat_sero_Japan_univariate$Abbott_ub <- NA

dat_sero_Japan_univariate$Roche_median <- NA
dat_sero_Japan_univariate$Roche_lb <- NA
dat_sero_Japan_univariate$Roche_ub <- NA

## Read in model to generate new quantities
model_weighted_Se <- stan_model(model_code=calc_wtd_Se_by_time_adjusting_for_severity_2groups)

## Generate severity-weighted sensitivities
weighted_Se_Abbott <- gqs(model_weighted_Se, data=list(prop=c(1-prop_hosp, prop_hosp)), draws=as.matrix(fit_Stan_Abbott))

weighted_Se_Roche <- gqs(model_weighted_Se, data=list(prop=c(1-prop_hosp, prop_hosp)), draws=as.matrix(fit_Stan_Roche))

## Get the inputs: time series
N <- dat_ts_Japan_final %>% nrow()
prop_infections <- dat_ts_Japan_final %>% select(prop_cases) %>% pull()

## Get the inputs: sensitivity
dat_Se_Abbott_Japan <- as.data.frame(weighted_Se_Abbott) %>% select(contains("wtd")) %>% select(1:N) %>% as.matrix()

dat_Se_Roche_Japan <- as.data.frame(weighted_Se_Roche) %>% select(contains("wtd")) %>% select(1:N) %>% as.matrix()

## Read in model to generate new quantities
model_weighted_Se_with_ts <- stan_model(model_code=calc_wtd_Se_by_time_adjusting_for_severity_and_ts)

## Estimate the parameters of Se (mean and sd), by assay
weighted_Se_Abbott_with_ts <- gqs(model_weighted_Se_with_ts, draws=dat_Se_Abbott_Japan, data=list(N=N, prop_infections=prop_infections))

weighted_Se_Roche_with_ts <- gqs(model_weighted_Se_with_ts, draws=dat_Se_Roche_Japan, data=list(N=N, prop_infections=prop_infections))

for(i in 1:nrow(dat_sero_Japan_univariate)) {
	
	## Fit the model - Abbott
	## Here: we use the binomial model because we know the numerators and denominators
	fit_Stan_Abbott <- stan(model_code=calc_seroprev_binomial,
	data=list(
		y_sample = dat_sero_Japan_univariate$Positive_Abbott[i] %>% as.integer(),
		n_sample = dat_sero_Japan_univariate$N_tested[i] %>% as.integer(),
		mu_Se = summary(weighted_Se_Abbott_with_ts)$summary["sensitivity_wtd_with_ts","mean"],
		sigma_Se = summary(weighted_Se_Abbott_with_ts)$summary["sensitivity_wtd_with_ts","sd"],
		Sp = Sp_Abbott),
	refresh=0)
	
	## Store
	 dat_sero_Japan_univariate$Abbott_median[i] <- summary(fit_Stan_Abbott)$summary["prev_adj","50%"]
	 dat_sero_Japan_univariate$Abbott_lb[i] <- summary(fit_Stan_Abbott)$summary["prev_adj","2.5%"]
	 dat_sero_Japan_univariate$Abbott_ub[i] <- summary(fit_Stan_Abbott)$summary["prev_adj","97.5%"]
	
	## Fit the model - Roche
	## Here: we use the binomial model because we know the numerators and denominators
	fit_Stan_Roche <- stan(model_code=calc_seroprev_binomial,
	data=list(
		y_sample = dat_sero_Japan_univariate$Positive_Roche[i] %>% as.integer(),
		n_sample = dat_sero_Japan_univariate$N_tested[i] %>% as.integer(),
		mu_Se = summary(weighted_Se_Roche_with_ts)$summary["sensitivity_wtd_with_ts","mean"],
		sigma_Se = summary(weighted_Se_Roche_with_ts)$summary["sensitivity_wtd_with_ts","sd"],
		Sp = Sp_Roche),
	refresh=0)
	
	## Store
	 dat_sero_Japan_univariate$Roche_median[i] <- summary(fit_Stan_Roche)$summary["prev_adj","50%"]
	 dat_sero_Japan_univariate$Roche_lb[i] <- summary(fit_Stan_Roche)$summary["prev_adj","2.5%"]
	 dat_sero_Japan_univariate$Roche_ub[i] <- summary(fit_Stan_Roche)$summary["prev_adj","97.5%"]
	
}

## Plot
windows()
dat_sero_Japan_univariate %>%
	pivot_longer(cols=Abbott_raw:Roche_ub) %>%
	separate(name, sep="_", into=c("Assay","metric")) %>%
	filter(metric=="raw") %>%
	ggplot() +
		geom_pointrange(data=dat_sero_Japan_univariate %>% mutate(Assay="Abbott") %>% mutate(metric="adjusted"), aes(x=Assay, y=Abbott_median, ymin=Abbott_lb, ymax=Abbott_ub, colour=metric)) +
		geom_pointrange(data=dat_sero_Japan_univariate %>% mutate(Assay="Roche") %>% mutate(metric="adjusted"), aes(x=Assay, y=Roche_median, ymin=Roche_lb, ymax=Roche_ub, colour=metric)) +
		geom_point(aes(x=Assay, y=value, colour=metric), size=2.5) +
		facet_wrap(.~Prefecture, scales="fixed", ncol=5) +
		scale_colour_manual(values=c("raw"="black", "adjusted"="red")) +
		theme_bw() +
		xlab("") +
		ylab("Seroprevalence")

## Prep to save
dat_sero_Japan_univariate_2groups <- dat_sero_Japan_univariate

## Save for future plotting etc
save(
	dat_sero_Japan_univariate_2groups,
	file="Results/Results_Japan_2groups_1assay.RData")
