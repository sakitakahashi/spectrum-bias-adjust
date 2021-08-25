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

FILENAME_VitrosIgG <- "Results/Stan_fit_VitrosIgG_2groups.RData"
load(FILENAME_VitrosIgG)
fit_Stan_VitrosIgG <- fit_Stan

## Read in model to generate new quantities
model_weighted_Se <- stan_model(model_code=calc_wtd_Se_by_time_adjusting_for_severity_2groups)

## Fix parameters
dat_severity_by_age <- read_excel("Data/Demography_data_US_severity_and_AS_with_N.xlsx")

## Get the estimates of specificity
Sp_Abbott <- 0.9963
Sp_Roche <- 0.9980
Sp_VitrosIgG <- 1.0

##############################################################################################
## Calculate severity-weighted & epi curve-weighted sensitivity for each assay, and adjust: ##
## 																						                                            	##
## US: using EpiNow2 & DEATHS															                                	##
##############################################################################################

## Read in merged EpiNow2 data
load("Results/EpiNow2_fit_US.RData")

## Add in census divisions
STORAGE_DEATHS_EpiNow2 %>%
	mutate(census_division = state) %>%
	mutate(census_division = ifelse(state %in% c("AK", "CA", "HI", "OR", "WA"), "Pacific", census_division)) %>%
	mutate(census_division = ifelse(state %in% c("AZ", "CO", "ID", "MT", "NV", "NM", "UT", "WY"), "Mountain", census_division)) %>%
	mutate(census_division = ifelse(state %in% c("IA", "KS", "MN", "MO", "NE", "ND", "SD"), "West North Central", census_division)) %>%
	mutate(census_division = ifelse(state %in% c("IL", "IN", "MI", "OH", "WI"), "East North Central", census_division)) %>%
	mutate(census_division = ifelse(state %in% c("AR", "LA", "OK", "TX"), "West South Central", census_division)) %>%
	mutate(census_division = ifelse(state %in% c("AL", "KY", "MS", "TN"), "East South Central", census_division)) %>%
	mutate(census_division = ifelse(state %in% c("NJ", "NY", "NYC", "PA"), "Middle Atlantic", census_division)) %>%
	mutate(census_division = ifelse(state %in% c("CT", "ME", "MA", "NH", "RI", "VT"), "New England", census_division)) %>%
	mutate(census_division = ifelse(state %in% c("DE", "FL", "GA", "MD", "NC", "SC", "VA", "DC", "WV", "PR"), "South Atlantic & Puerto Rico", census_division)) %>%
	filter(census_division %in% c("Pacific", "Mountain", "West North Central", "East North Central", "West South Central", "East South Central", "Middle Atlantic", "New England", "South Atlantic & Puerto Rico")) %>%
	rename(Time = summarised.date) %>%
	rename(onsets = summarised.median) %>%
	droplevels() %>%
	select(Time, onsets, state, census_division) -> data_frame_all_states_NYT_DEATHS

################################################################################

## FUNCTION
calc_wtd_Se_by_time_adjusting_for_severity_and_ts_2_assays_propAssay <-'
data {
	
	int N;
	real prop_infections[N];
	real prop_Abbott;
	
}
parameters {
	
	real sensitivity_wtd_Abbott[N];
	real sensitivity_wtd_Roche[N];
	
}
generated quantities {
	
	real sensitivity_wtd_Abbott_with_ts;
	real sensitivity_wtd_Roche_with_ts;
	
	real sensitivity_wtd_overall_with_ts;
	
	sensitivity_wtd_Abbott_with_ts = 0.0;
	sensitivity_wtd_Roche_with_ts = 0.0;
	
	for(i in 1:N) {
		
		sensitivity_wtd_Abbott_with_ts += sensitivity_wtd_Abbott[i] * prop_infections[N+1-i];
		sensitivity_wtd_Roche_with_ts += sensitivity_wtd_Roche[i] * prop_infections[N+1-i];
		
	}
	
	sensitivity_wtd_overall_with_ts = (prop_Abbott*sensitivity_wtd_Abbott_with_ts) + (1.0-prop_Abbott)*sensitivity_wtd_Roche_with_ts;
	
}
'

## Read in model to generate new quantities (2 assays)
model_weighted_Se_with_ts <- stan_model(model_code=calc_wtd_Se_by_time_adjusting_for_severity_and_ts_2_assays_propAssay)

## Pull the severity by age
## Assuming doesn't vary by sex
## Weight it by age and by census division to get a single value
dat_severity_by_age %>%
	group_by(census_division) %>%
	summarise(
		prop_hosp = weighted.mean(prob_hosp, N),
		prop_AS = weighted.mean(prob_AS, N),
		N = sum(N)) %>%
	summarise(prop_hosp_wtd = weighted.mean(prop_hosp, N)) %>%
	select(prop_hosp_wtd) %>%
	pull() -> prop_hosp_wtd

#####################################
## LOOP OVER DATES AND PROP_ABBOTT ##
#####################################

dat_track_US_assay_time <- expand.grid(date=seq(ymd("2020-03-31"), ymd("2021-03-15"), by="2 weeks") - weeks(3), prop_Abbott=seq(0, 1, by=0.05), census_division=unique(data_frame_all_states_NYT_DEATHS$census_division))
dat_track_US_assay_time$wtd_Se <- NA

dim(dat_track_US_assay_time)

## ## Generate severity-weighted sensitivities
## Leaving out Vitros
weighted_Se_Abbott <- gqs(model_weighted_Se, data=list(prop=c(1-prop_hosp_wtd, prop_hosp_wtd)), draws=as.matrix(fit_Stan_Abbott))

weighted_Se_Roche <- gqs(model_weighted_Se, data=list(prop=c(1-prop_hosp_wtd, prop_hosp_wtd)), draws=as.matrix(fit_Stan_Roche))

for(xx in 1:nrow(dat_track_US_assay_time)) {
	
	## Which dates?
	data_frame_all_states_NYT_DEATHS %>%
		filter(Time <= dat_track_US_assay_time$date[xx]) %>%
		mutate(days_before_survey = as.double(dat_track_US_assay_time$date[xx] - Time)) %>%
		ungroup() %>%
		filter(census_division==dat_track_US_assay_time$census_division[xx]) %>%
		group_by(census_division, Time) %>%
		summarise(onsets = sum(onsets)) %>%
		ungroup() %>%
		mutate(onsets_prop = onsets / sum(onsets)) %>%
		pull() -> PROP_INFECTIONS
	
	N <- length(PROP_INFECTIONS)
	
	## Get the inputs: sensitivity
	dat_Se_Abbott_US <- as.data.frame(weighted_Se_Abbott) %>% select(contains("wtd")) %>% select(1:N)
	dat_Se_Roche_US <- as.data.frame(weighted_Se_Roche) %>% select(contains("wtd")) %>% select(1:N)
	
	names(dat_Se_Abbott_US) <- dat_Se_Abbott_US %>% names() %>% str_match("(sensitivity_wtd)(\\[\\d+\\])") %>% {paste0(.[,2],"_Abbott", .[,3])}
	names(dat_Se_Roche_US) <- dat_Se_Roche_US %>% names() %>% str_match("(sensitivity_wtd)(\\[\\d+\\])") %>% {paste0(.[,2],"_Roche", .[,3])}
	
	## Collate the data from 2 assays
	dat_US_kinetics <- cbind(as.matrix(dat_Se_Abbott_US), as.matrix(dat_Se_Roche_US))
	
	## Estimate the parameters of Se (mean and sd), by assay
	weighted_Se_with_ts <- gqs(model_weighted_Se_with_ts, draws=dat_US_kinetics, data=list(N=N, prop_infections=PROP_INFECTIONS, prop_Abbott=dat_track_US_assay_time$prop_Abbott[xx]))
	
	## Save
	dat_track_US_assay_time$wtd_Se[xx] <- summary(weighted_Se_with_ts)$summary["sensitivity_wtd_overall_with_ts","50%"]
	
	print(xx)
	
}

## Save for future plotting etc
save(data_frame_all_states_NYT_DEATHS, dat_track_US_assay_time, file="Results/Results_US_2groups_loopOverAssayAndTime.RData")

################################################################################

load("Results/Results_US_2groups_loopOverAssayAndTime.RData")

data_frame_all_states_NYT_DEATHS %>%
  filter(Time <= ymd("2021-02-20")) %>%
  ggplot(aes(x=Time, y=onsets, group=state)) +
  geom_line(size=0.1) +
  scale_x_date(date_breaks="1 month") +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(.~state, scales="free_y") +
  ggtitle("EpiNow2 (deaths)")

data_frame_all_states_NYT_DEATHS %>%
  filter(Time <= ymd("2021-02-20")) %>%
  group_by(census_division, Time) %>%
  summarise(onsets = sum(onsets)) %>%
  ungroup() %>%
  mutate(onsets_prop = onsets / sum(onsets)) %>%
  rename(date = Time) %>%
  as.data.frame() %>%
  mutate(census_division = factor(as.character(census_division))) %>%
  group_by(census_division) %>%
  mutate(onsets_prop_scaled = onsets_prop / max(onsets_prop)) %>%
  ungroup() -> ts_to_overlay

dat_track_US_assay_time %>%
  mutate(census_division = factor(as.character(census_division))) %>%
  ggplot(aes(x=date, y=prop_Abbott)) +
  geom_tile(aes(fill=wtd_Se)) +
  geom_contour(aes(z=wtd_Se), colour="darkgrey", binwidth=0.05) +
  geom_line(data=ts_to_overlay, aes(x=date, y=onsets_prop_scaled), colour="black", size=1.1) +
  theme_bw() +
  ylab("Proportion of tests on Abbott ARCHITECT\n(assuming the rest are on Roche Elecsys)") +
  xlab("Date of serosurvey") +
  scale_x_date(date_breaks="2 months") +
  scale_y_continuous(breaks=scales::pretty_breaks(6)) +
  scale_fill_distiller("Weighted\nsensitivity", palette="RdYlBu", direction=1) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  facet_wrap(.~census_division, ncol=3, scales="fixed") +
  ggtitle("United States (2 severity groups)")
