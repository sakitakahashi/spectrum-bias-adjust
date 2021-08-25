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

## THIS VERSION GENERATES THE AGE-ADJUSTED ESTIMATES (as opposed to stratified)

source("Code/functions_for_gqs.R")

########################
## Load kinetics data ##
########################

## Load each assay
FILENAME_Abbott <- "Results/Stan_fit_Abbott_3groups.RData"
load(FILENAME_Abbott)
fit_Stan_Abbott <- fit_Stan

FILENAME_Roche <- "Results/Stan_fit_Roche_3groups.RData"
load(FILENAME_Roche)
fit_Stan_Roche <- fit_Stan

FILENAME_VitrosIgG <- "Results/Stan_fit_VitrosIgG_3groups.RData"
load(FILENAME_VitrosIgG)
fit_Stan_VitrosIgG <- fit_Stan

## Read in model to generate new quantities
## All this is done below
model_weighted_Se <- stan_model(model_code=calc_wtd_Se_by_time_adjusting_for_severity)

## Fix parameters
dat_severity_by_age <- read_excel("Data/Demography_data_US_severity_and_AS.xlsx")

## Get the estimates of specificity
Sp_Abbott <- 0.9963
Sp_Roche <- 0.9980
Sp_VitrosIgG <- 1.0

##############################################################################################
## Calculate severity-weighted & epi curve-weighted sensitivity for each assay, and adjust: ##
## 																						                                            	##
## US: using EpiNow2 & deaths														                                		##
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

ex <- interval(ymd("2020-07-27"), ymd("2020-08-13"))
int_start(ex) + (int_end(ex) - int_start(ex))/2

ex <- interval(ymd("2020-08-10"), ymd("2020-08-27"))
int_start(ex) + (int_end(ex) - int_start(ex))/2

ex <- interval(ymd("2020-08-24"), ymd("2020-09-10"))
int_start(ex) + (int_end(ex) - int_start(ex))/2

ex <- interval(ymd("2020-09-07"), ymd("2020-09-24"))
int_start(ex) + (int_end(ex) - int_start(ex))/2

## Remove 3 weeks
SURVEY_DATE_US_1 <- ymd("2020-08-04")-weeks(3)
SURVEY_DATE_US_2 <- ymd("2020-08-18")-weeks(3)
SURVEY_DATE_US_3 <- ymd("2020-09-01")-weeks(3)
SURVEY_DATE_US_4 <- ymd("2020-09-15")-weeks(3)

## Which dates?
data_frame_all_states_NYT_DEATHS %>%
	filter(Time <= SURVEY_DATE_US_1) %>%
	mutate(days_before_survey_1 = as.double(SURVEY_DATE_US_1 - Time)) %>%
	ungroup() -> dat_ts_US_1_final

data_frame_all_states_NYT_DEATHS %>%
	filter(Time <= SURVEY_DATE_US_2) %>%
	mutate(days_before_survey_2 = as.double(SURVEY_DATE_US_2 - Time)) %>%
	ungroup() -> dat_ts_US_2_final

data_frame_all_states_NYT_DEATHS %>%
	filter(Time <= SURVEY_DATE_US_3) %>%
	mutate(days_before_survey_3 = as.double(SURVEY_DATE_US_3 - Time)) %>%
	ungroup() -> dat_ts_US_3_final

data_frame_all_states_NYT_DEATHS %>%
	filter(Time <= SURVEY_DATE_US_4) %>%
	mutate(days_before_survey_4 = as.double(SURVEY_DATE_US_4 - Time)) %>%
	ungroup() -> dat_ts_US_4_final

############################
## ACTUAL ADJUSTMENT PART ##
############################

## Get number of tests per date range, per state, per round
## Not publicly available
dat_US_state_counts <- read_csv(...)

dat_US_state_counts %>%
	rename(state = Site) %>%
	rename(Date = `Date Range of Specimen Collection`) %>%
	rename(n_Age_0_17 = `n [0-17 Years Prevalence]`) %>%
	rename(n_Age_18_49 = `n [18-49 Years Prevalence]`) %>%
	rename(n_Age_50_64 = `n [50-64 Years Prevalence]`) %>%
	rename(n_Age_65_plus = `n [65+ Years Prevalence]`) %>%
	rename(n_Female = `n [Female Prevalence]`) %>%
	rename(n_Male = `n [Male Prevalence]`) %>%
	rename(n_Total = `n [Cumulative Prevalence]`) %>%
	select(state, Round, Date, n_Age_0_17, n_Age_18_49, n_Age_50_64, n_Age_65_plus, n_Female, n_Male, n_Total) %>%
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
	mutate(census_division = factor(census_division, levels=c("Pacific", "Mountain", "West North Central", "East North Central", "West South Central", "East South Central", "Middle Atlantic", "New England", "South Atlantic & Puerto Rico"))) %>%
	arrange(Round, census_division) -> dat_US_state_counts

## Get proportion by group
dat_US_state_counts %>%
	group_by(census_division, Round) %>%
	mutate(prop_Age_0_17 = n_Age_0_17 / sum(n_Age_0_17)) %>%
	mutate(prop_Age_18_49 = n_Age_18_49 / sum(n_Age_18_49)) %>%
	mutate(prop_Age_50_64 = n_Age_50_64 / sum(n_Age_50_64)) %>%
	mutate(prop_Age_65_plus = n_Age_65_plus / sum(n_Age_65_plus)) %>%
	ungroup() -> dat_US_state_counts

## Generate weighted time series
## Didn't do this by sex because we don't have all the pairs (sex-age)
dat_ts_US_1_final %>%
	left_join(dat_US_state_counts %>% filter(Round==1)) %>%
	group_by(Time, days_before_survey_1, census_division) %>%
	summarise(
		onsets_mean_Age_0_17 = weighted.mean(onsets, prop_Age_0_17),
		onsets_mean_Age_18_49 = weighted.mean(onsets, prop_Age_18_49),
		onsets_mean_Age_50_64 = weighted.mean(onsets, prop_Age_50_64),
		onsets_mean_Age_65_plus = weighted.mean(onsets, prop_Age_65_plus)) %>%
	## NEED TO DO THIS BECAUSE SOME COUNTS ARE LOW
	mutate(
		onsets_mean_Age_0_17 = ifelse(is.na(onsets_mean_Age_0_17),0,onsets_mean_Age_0_17),
		onsets_mean_Age_18_49 = ifelse(is.na(onsets_mean_Age_18_49),0,onsets_mean_Age_18_49),
		onsets_mean_Age_50_64 = ifelse(is.na(onsets_mean_Age_50_64),0,onsets_mean_Age_50_64),
		onsets_mean_Age_65_plus = ifelse(is.na(onsets_mean_Age_65_plus),0,onsets_mean_Age_65_plus)) %>%
	group_by(census_division) %>%
	mutate(
		onsets_prop_Age_0_17 = onsets_mean_Age_0_17 / sum(onsets_mean_Age_0_17),
		onsets_prop_Age_18_49 = onsets_mean_Age_18_49 / sum(onsets_mean_Age_18_49),
		onsets_prop_Age_50_64 = onsets_mean_Age_50_64 / sum(onsets_mean_Age_50_64),
		onsets_prop_Age_65_plus = onsets_mean_Age_65_plus / sum(onsets_mean_Age_65_plus)) %>%
	ungroup() -> dat_ts_US_1_final

dat_ts_US_2_final %>%
	left_join(dat_US_state_counts %>% filter(Round==2)) %>%
	group_by(Time, days_before_survey_2, census_division) %>%
	summarise(
		onsets_mean_Age_0_17 = weighted.mean(onsets, prop_Age_0_17),
		onsets_mean_Age_18_49 = weighted.mean(onsets, prop_Age_18_49),
		onsets_mean_Age_50_64 = weighted.mean(onsets, prop_Age_50_64),
		onsets_mean_Age_65_plus = weighted.mean(onsets, prop_Age_65_plus)) %>%
	## NEED TO DO THIS BECAUSE SOME COUNTS ARE LOW
	mutate(
		onsets_mean_Age_0_17 = ifelse(is.na(onsets_mean_Age_0_17),0,onsets_mean_Age_0_17),
		onsets_mean_Age_18_49 = ifelse(is.na(onsets_mean_Age_18_49),0,onsets_mean_Age_18_49),
		onsets_mean_Age_50_64 = ifelse(is.na(onsets_mean_Age_50_64),0,onsets_mean_Age_50_64),
		onsets_mean_Age_65_plus = ifelse(is.na(onsets_mean_Age_65_plus),0,onsets_mean_Age_65_plus)) %>%
	group_by(census_division) %>%
	mutate(
		onsets_prop_Age_0_17 = onsets_mean_Age_0_17 / sum(onsets_mean_Age_0_17),
		onsets_prop_Age_18_49 = onsets_mean_Age_18_49 / sum(onsets_mean_Age_18_49),
		onsets_prop_Age_50_64 = onsets_mean_Age_50_64 / sum(onsets_mean_Age_50_64),
		onsets_prop_Age_65_plus = onsets_mean_Age_65_plus / sum(onsets_mean_Age_65_plus)) %>%
	ungroup() -> dat_ts_US_2_final

dat_ts_US_3_final %>%
	left_join(dat_US_state_counts %>% filter(Round==3)) %>%
	group_by(Time, days_before_survey_3, census_division) %>%
	summarise(
		onsets_mean_Age_0_17 = weighted.mean(onsets, prop_Age_0_17),
		onsets_mean_Age_18_49 = weighted.mean(onsets, prop_Age_18_49),
		onsets_mean_Age_50_64 = weighted.mean(onsets, prop_Age_50_64),
		onsets_mean_Age_65_plus = weighted.mean(onsets, prop_Age_65_plus)) %>%
	## NEED TO DO THIS BECAUSE SOME COUNTS ARE LOW
	mutate(
		onsets_mean_Age_0_17 = ifelse(is.na(onsets_mean_Age_0_17),0,onsets_mean_Age_0_17),
		onsets_mean_Age_18_49 = ifelse(is.na(onsets_mean_Age_18_49),0,onsets_mean_Age_18_49),
		onsets_mean_Age_50_64 = ifelse(is.na(onsets_mean_Age_50_64),0,onsets_mean_Age_50_64),
		onsets_mean_Age_65_plus = ifelse(is.na(onsets_mean_Age_65_plus),0,onsets_mean_Age_65_plus)) %>%
	group_by(census_division) %>%
	mutate(
		onsets_prop_Age_0_17 = onsets_mean_Age_0_17 / sum(onsets_mean_Age_0_17),
		onsets_prop_Age_18_49 = onsets_mean_Age_18_49 / sum(onsets_mean_Age_18_49),
		onsets_prop_Age_50_64 = onsets_mean_Age_50_64 / sum(onsets_mean_Age_50_64),
		onsets_prop_Age_65_plus = onsets_mean_Age_65_plus / sum(onsets_mean_Age_65_plus)) %>%
	ungroup() -> dat_ts_US_3_final

dat_ts_US_4_final %>%
	left_join(dat_US_state_counts %>% filter(Round==4)) %>%
	group_by(Time, days_before_survey_4, census_division) %>%
	summarise(
		onsets_mean_Age_0_17 = weighted.mean(onsets, prop_Age_0_17),
		onsets_mean_Age_18_49 = weighted.mean(onsets, prop_Age_18_49),
		onsets_mean_Age_50_64 = weighted.mean(onsets, prop_Age_50_64),
		onsets_mean_Age_65_plus = weighted.mean(onsets, prop_Age_65_plus)) %>%
	## NEED TO DO THIS BECAUSE SOME COUNTS ARE LOW
	mutate(
		onsets_mean_Age_0_17 = ifelse(is.na(onsets_mean_Age_0_17),0,onsets_mean_Age_0_17),
		onsets_mean_Age_18_49 = ifelse(is.na(onsets_mean_Age_18_49),0,onsets_mean_Age_18_49),
		onsets_mean_Age_50_64 = ifelse(is.na(onsets_mean_Age_50_64),0,onsets_mean_Age_50_64),
		onsets_mean_Age_65_plus = ifelse(is.na(onsets_mean_Age_65_plus),0,onsets_mean_Age_65_plus)) %>%
	group_by(census_division) %>%
	mutate(
		onsets_prop_Age_0_17 = onsets_mean_Age_0_17 / sum(onsets_mean_Age_0_17),
		onsets_prop_Age_18_49 = onsets_mean_Age_18_49 / sum(onsets_mean_Age_18_49),
		onsets_prop_Age_50_64 = onsets_mean_Age_50_64 / sum(onsets_mean_Age_50_64),
		onsets_prop_Age_65_plus = onsets_mean_Age_65_plus / sum(onsets_mean_Age_65_plus)) %>%
	ungroup() -> dat_ts_US_4_final

#############################
## Bring in the serosurvey ##
#############################

## Bring in serosurvey results
## Not publicly available
dat_sero_US <- read_csv(...)

dat_sero_US %>%
	mutate(Round = Testing_Period) %>%
	mutate(Testing_Period = factor(Testing_Period, levels=1:4, labels=c("7/27-8/13", "8/10-8/27", "8/24-9/10", "9/7-9/24"))) %>%
	mutate(Sex = factor(Sex, levels=1:2, labels=c("M", "F"))) %>%
	mutate(Age_Category = factor(Age_Category, levels=1:4, labels=c("0-17", "18-49", "50-64", "65+"))) %>%
	mutate(Assay = factor(Assay, levels=2:4, labels=c("Roche Elecsys", "Abbott ARCHITECT", "Ortho VITROS IgG"))) %>%
	mutate(Test_Result = factor(Test_Result, levels=0:1, labels=c("Negative", "Positive"))) %>%
	mutate(Census_Division = factor(Census_Division, levels=1:9, labels=c("Pacific", "Mountain", "West North Central", "East North Central", "West South Central", "East South Central", "Middle Atlantic", "New England", "South Atlantic & Puerto Rico"))) %>%
	as_tibble() %>%
	mutate_if(is.factor, forcats::fct_explicit_na) -> dat_sero_US

## Remove the NA age and sex classes, for now
dat_sero_US %>%
	filter(Age_Category!="(Missing)") %>%
	filter(Sex!="(Missing)") %>%
	droplevels() -> dat_sero_US

## Aggregate
dat_sero_US %>%
	group_by(Round, Testing_Period, Sex, Census_Division, Age_Category, Assay) %>%
	summarise(
		n_positive = sum(Test_Result=="Positive"),
		n_tested = n(),
		seroprevalence = n_positive/ n_tested) -> dat_sero_US

dim(dat_sero_US)

## For estimation: have 1 row for each assay
dat_sero_US %>%
	complete(Assay, nesting(Round, Testing_Period, Census_Division, Age_Category, Sex), fill=list(n_positive=0, n_tested=0, seroprevalence=0)) %>%
	group_by(Round, Testing_Period, Census_Division, Age_Category, Sex) %>%
	mutate(prop_tests = n_tested / sum(n_tested)) %>%
	ungroup() %>%
	as.data.frame() -> dat_sero_US

## Should be: 2*3*4*9*4 = 864
dim(dat_sero_US)

## For tracking: have 1 row for each point (across assays)
## Should be 4*9=36
dat_sero_US %>%
	group_by(Round, Testing_Period, Census_Division) %>%
	summarise(
		n_positive_raw = sum(n_positive),
		n_tested_raw = sum(n_tested),
		seroprevalence_raw = n_positive_raw / n_tested_raw) %>%
	ungroup() -> dat_sero_US_track

dat_sero_US_track$adj_sero_median <- NA
dat_sero_US_track$adj_sero_lb <- NA
dat_sero_US_track$adj_sero_ub <- NA

## Bring in demographic data
dat_pop_by_div <- read_csv("Data/Demography_data_US_age_4binstotal.csv")

dat_pop_by_div %>%
	select(Age, contains("!!Female!!Estimate"), contains("!!Male!!Estimate")) %>%
	as.data.frame() -> dat_pop_by_div

dat_pop_by_div$Age_Category <- c("0-17","18-49","18-49","18-49","18-49","18-49","18-49","50-64","50-64","50-64","65+")

dat_pop_by_div %>%
	pivot_longer(cols=contains("!!"), names_to="census_division_sex", values_to="N") %>%
	separate(census_division_sex, c("Census_Division","Sex","to_rm"), sep='!!') %>%
	select(-to_rm) %>%
	mutate(Sex = as.character(Sex)) %>%
	mutate(Sex = str_sub(Sex, 1,1)) %>%
	mutate(Census_Division = str_extract(Census_Division, ".+?(?= Division)")) %>%
	mutate(Census_Division = factor(Census_Division, levels=c("Pacific", "Mountain", "West North Central", "East North Central", "West South Central", "East South Central", "Middle Atlantic", "New England", "South Atlantic"), labels=c("Pacific", "Mountain", "West North Central", "East North Central", "West South Central", "East South Central", "Middle Atlantic", "New England", "South Atlantic & Puerto Rico"))) %>%
	group_by(Age_Category, Sex, Census_Division) %>%
	summarise(pop_of_census_division = sum(N)) %>%
	ungroup() %>%
	arrange(Census_Division, Age_Category, Sex) -> dat_pop_by_div_clean

dat_sero_US %>%
	left_join(dat_pop_by_div_clean, by=c("Census_Division","Age_Category","Sex")) %>%
	arrange(Round, Census_Division, Age_Category, Sex, Assay) -> dat_sero_US

## Read in model to generate new quantities (3 assays)
model_weighted_Se_with_ts <- stan_model(model_code=calc_wtd_Se_by_time_adjusting_for_severity_and_ts_3_assays)

for(i in 1:nrow(dat_sero_US_track)) {
	
	## Pull the survey meta-data
	WHICH_ROUND <- dat_sero_US_track$Round[i]
	WHICH_CENSUS_DIVISION <- dat_sero_US_track$Census_Division[i] %>% as.character()
	
	## Pull the severity by age
	## Assuming doesn't vary by sex
	dat_severity_by_age %>%
		filter(census_division==WHICH_CENSUS_DIVISION) %>%
		select(prob_hosp) %>%
		pull() -> prop_hosp
	
	## Pull the prop AS by age
	## Assuming doesn't vary by sex
	dat_severity_by_age %>%
		filter(census_division==WHICH_CENSUS_DIVISION) %>%
		select(prob_AS) %>%
		pull() -> prop_AS
	
	## Generate severity-weighted sensitivities
	weighted_Se_Abbott_1_0_to_17 <- gqs(model_weighted_Se, data=list(prop=c(prop_AS[1], 1-prop_AS[1]-prop_hosp[1], prop_hosp[1])), draws=as.matrix(fit_Stan_Abbott))
	weighted_Se_Abbott_2_18_to_49 <- gqs(model_weighted_Se, data=list(prop=c(prop_AS[2], 1-prop_AS[2]-prop_hosp[2], prop_hosp[2])), draws=as.matrix(fit_Stan_Abbott))
	weighted_Se_Abbott_3_50_to_64 <- gqs(model_weighted_Se, data=list(prop=c(prop_AS[3], 1-prop_AS[3]-prop_hosp[3], prop_hosp[3])), draws=as.matrix(fit_Stan_Abbott))
	weighted_Se_Abbott_4_65_plus <- gqs(model_weighted_Se, data=list(prop=c(prop_AS[4], 1-prop_AS[4]-prop_hosp[4], prop_hosp[4])), draws=as.matrix(fit_Stan_Abbott))

	weighted_Se_Roche_1_0_to_17 <- gqs(model_weighted_Se, data=list(prop=c(prop_AS[1], 1-prop_AS[1]-prop_hosp[1], prop_hosp[1])), draws=as.matrix(fit_Stan_Roche))
	weighted_Se_Roche_2_18_to_49 <- gqs(model_weighted_Se, data=list(prop=c(prop_AS[2], 1-prop_AS[2]-prop_hosp[2], prop_hosp[2])), draws=as.matrix(fit_Stan_Roche))
	weighted_Se_Roche_3_50_to_64 <- gqs(model_weighted_Se, data=list(prop=c(prop_AS[3], 1-prop_AS[3]-prop_hosp[3], prop_hosp[3])), draws=as.matrix(fit_Stan_Roche))
	weighted_Se_Roche_4_65_plus <- gqs(model_weighted_Se, data=list(prop=c(prop_AS[4], 1-prop_AS[4]-prop_hosp[4], prop_hosp[4])), draws=as.matrix(fit_Stan_Roche))
	
	weighted_Se_VitrosIgG_1_0_to_17 <- gqs(model_weighted_Se, data=list(prop=c(prop_AS[1], 1-prop_AS[1]-prop_hosp[1], prop_hosp[1])), draws=as.matrix(fit_Stan_VitrosIgG))
	weighted_Se_VitrosIgG_2_18_to_49 <- gqs(model_weighted_Se, data=list(prop=c(prop_AS[2], 1-prop_AS[2]-prop_hosp[2], prop_hosp[2])), draws=as.matrix(fit_Stan_VitrosIgG))
	weighted_Se_VitrosIgG_3_50_to_64 <- gqs(model_weighted_Se, data=list(prop=c(prop_AS[3], 1-prop_AS[3]-prop_hosp[3], prop_hosp[3])), draws=as.matrix(fit_Stan_VitrosIgG))
	weighted_Se_VitrosIgG_4_65_plus <- gqs(model_weighted_Se, data=list(prop=c(prop_AS[4], 1-prop_AS[4]-prop_hosp[4], prop_hosp[4])), draws=as.matrix(fit_Stan_VitrosIgG))
	
	rm(prop_hosp, prop_AS)
	
	## Pull the serosurvey results
	## Order: Abbott, Ortho, Roche
	# dat_sero_US %>%
		# filter(Round==WHICH_ROUND & Census_Division==WHICH_CENSUS_DIVISION & Age_Category==WHICH_AGE_CATEGORY) %>%
		# droplevels() %>%
		# mutate(Assay = as.character(Assay)) %>%
		# arrange(Assay) %>%
		# select(seroprevalence) %>%
		# pull() -> SEROPREVALENCE_BY_ASSAY
	
	###### Get raw data
	## Assays: Abbott ARCHITECT, Ortho VITROS IgG, Roche Elecsys
	
	## [01] F, 0-17
	## [02] M, 0-17
	## [03] F, 18-49
	## [04] M, 18-49
	## [05] F, 50-64
	## [06] M, 50-64
	## [07] F, 65+
	## [08] M, 65+
	
	## Separated by assay
	
	## Abbott
	dat_sero_US %>%
		filter(Round==WHICH_ROUND & Census_Division==WHICH_CENSUS_DIVISION) %>%
		droplevels() %>%
		mutate(Assay = as.character(Assay)) %>%
		arrange(Age_Category, Sex, Assay) %>%
		filter(Assay == "Abbott ARCHITECT") -> tmp_Abbott

	tmp_Abbott %>%
		select(n_positive) %>%
		pull() -> NUMERATOR_BY_ASSAY_Abbott

	tmp_Abbott %>%
		select(n_tested) %>%
		pull() -> DENOMINATOR_BY_ASSAY_Abbott

	dat_sero_US %>%
		filter(Round==WHICH_ROUND & Census_Division==WHICH_CENSUS_DIVISION) %>%
		droplevels() %>%
		mutate(Assay = as.character(Assay)) %>%
		arrange(Age_Category, Sex, Assay) %>%
		filter(Assay == "Ortho VITROS IgG") -> tmp_Ortho
	
	## Ortho
	tmp_Ortho %>%
		select(n_positive) %>%
		pull() -> NUMERATOR_BY_ASSAY_Ortho
	
	tmp_Ortho %>%
		select(n_tested) %>%
		pull() -> DENOMINATOR_BY_ASSAY_Ortho
	
	## Roche
	dat_sero_US %>%
		filter(Round==WHICH_ROUND & Census_Division==WHICH_CENSUS_DIVISION) %>%
		droplevels() %>%
		mutate(Assay = as.character(Assay)) %>%
		arrange(Age_Category, Sex, Assay) %>%
		filter(Assay == "Roche Elecsys") -> tmp_Roche
	
	tmp_Roche %>%
		select(n_positive) %>%
		pull() -> NUMERATOR_BY_ASSAY_Roche
	
	tmp_Roche %>%
		select(n_tested) %>%
		pull() -> DENOMINATOR_BY_ASSAY_Roche
	
	# dat_sero_US %>%
		# filter(Round==WHICH_ROUND & Census_Division==WHICH_CENSUS_DIVISION) %>%
		# droplevels() %>%
		# mutate(Assay = as.character(Assay)) %>%
		# arrange(Age_Category, Sex, Assay) %>%
		# select(prop_tests) %>%
		# pull() -> PROP_TEST_BY_ASSAY
	
	## Pull the correct time series
	if(WHICH_ROUND==1) DAT_FOR_NOW <- dat_ts_US_1_final
	if(WHICH_ROUND==2) DAT_FOR_NOW <- dat_ts_US_2_final
	if(WHICH_ROUND==3) DAT_FOR_NOW <- dat_ts_US_3_final
	if(WHICH_ROUND==4) DAT_FOR_NOW <- dat_ts_US_4_final
	
	DAT_FOR_NOW %>%
		filter(census_division==WHICH_CENSUS_DIVISION) %>%
		select(onsets_prop_Age_0_17) %>%
		pull() -> PROP_INFECTIONS_age_1_0_to_17
	
	DAT_FOR_NOW %>%
		filter(census_division==WHICH_CENSUS_DIVISION) %>%
		select(onsets_prop_Age_18_49) %>%
		pull() -> PROP_INFECTIONS_age_2_18_to_49
	
	DAT_FOR_NOW %>%
		filter(census_division==WHICH_CENSUS_DIVISION) %>%
		select(onsets_prop_Age_50_64) %>%
		pull() -> PROP_INFECTIONS_age_3_50_to_64

	DAT_FOR_NOW %>%
		filter(census_division==WHICH_CENSUS_DIVISION) %>%
		select(onsets_prop_Age_65_plus) %>%
		pull() -> PROP_INFECTIONS_age_4_65_plus
	
	N <- length(PROP_INFECTIONS_age_1_0_to_17)
	
	## Get the inputs: sensitivity
	## THIS DEPENDS ON AGE
	dat_Se_Abbott_US_age_1 <- as.data.frame(weighted_Se_Abbott_1_0_to_17) %>% select(contains("wtd")) %>% select(1:N)
	dat_Se_Abbott_US_age_2 <- as.data.frame(weighted_Se_Abbott_2_18_to_49) %>% select(contains("wtd")) %>% select(1:N)
	dat_Se_Abbott_US_age_3 <- as.data.frame(weighted_Se_Abbott_3_50_to_64) %>% select(contains("wtd")) %>% select(1:N)
	dat_Se_Abbott_US_age_4 <- as.data.frame(weighted_Se_Abbott_4_65_plus) %>% select(contains("wtd")) %>% select(1:N)
	
	dat_Se_VitrosIgG_US_age_1 <- as.data.frame(weighted_Se_VitrosIgG_1_0_to_17) %>% select(contains("wtd")) %>% select(1:N)
	dat_Se_VitrosIgG_US_age_2 <- as.data.frame(weighted_Se_VitrosIgG_2_18_to_49) %>% select(contains("wtd")) %>% select(1:N)
	dat_Se_VitrosIgG_US_age_3 <- as.data.frame(weighted_Se_VitrosIgG_3_50_to_64) %>% select(contains("wtd")) %>% select(1:N)
	dat_Se_VitrosIgG_US_age_4 <- as.data.frame(weighted_Se_VitrosIgG_4_65_plus) %>% select(contains("wtd")) %>% select(1:N)
	
	dat_Se_Roche_US_age_1 <- as.data.frame(weighted_Se_Roche_1_0_to_17) %>% select(contains("wtd")) %>% select(1:N)
	dat_Se_Roche_US_age_2 <- as.data.frame(weighted_Se_Roche_2_18_to_49) %>% select(contains("wtd")) %>% select(1:N)
	dat_Se_Roche_US_age_3 <- as.data.frame(weighted_Se_Roche_3_50_to_64) %>% select(contains("wtd")) %>% select(1:N)
	dat_Se_Roche_US_age_4 <- as.data.frame(weighted_Se_Roche_4_65_plus) %>% select(contains("wtd")) %>% select(1:N)
	
	names(dat_Se_Abbott_US_age_1) <- dat_Se_Abbott_US_age_1 %>% names() %>% str_match("(sensitivity_wtd)(\\[\\d+\\])") %>% {paste0(.[,2],"_Abbott", .[,3])}
	names(dat_Se_Abbott_US_age_2) <- dat_Se_Abbott_US_age_2 %>% names() %>% str_match("(sensitivity_wtd)(\\[\\d+\\])") %>% {paste0(.[,2],"_Abbott", .[,3])}
	names(dat_Se_Abbott_US_age_3) <- dat_Se_Abbott_US_age_3 %>% names() %>% str_match("(sensitivity_wtd)(\\[\\d+\\])") %>% {paste0(.[,2],"_Abbott", .[,3])}
	names(dat_Se_Abbott_US_age_4) <- dat_Se_Abbott_US_age_4 %>% names() %>% str_match("(sensitivity_wtd)(\\[\\d+\\])") %>% {paste0(.[,2],"_Abbott", .[,3])}
	
	names(dat_Se_VitrosIgG_US_age_1) <- dat_Se_VitrosIgG_US_age_1 %>% names() %>% str_match("(sensitivity_wtd)(\\[\\d+\\])") %>% {paste0(.[,2],"_VitrosIgG", .[,3])}
	names(dat_Se_VitrosIgG_US_age_2) <- dat_Se_VitrosIgG_US_age_2 %>% names() %>% str_match("(sensitivity_wtd)(\\[\\d+\\])") %>% {paste0(.[,2],"_VitrosIgG", .[,3])}
	names(dat_Se_VitrosIgG_US_age_3) <- dat_Se_VitrosIgG_US_age_3 %>% names() %>% str_match("(sensitivity_wtd)(\\[\\d+\\])") %>% {paste0(.[,2],"_VitrosIgG", .[,3])}
	names(dat_Se_VitrosIgG_US_age_4) <- dat_Se_VitrosIgG_US_age_4 %>% names() %>% str_match("(sensitivity_wtd)(\\[\\d+\\])") %>% {paste0(.[,2],"_VitrosIgG", .[,3])}
	
	names(dat_Se_Roche_US_age_1) <- dat_Se_Roche_US_age_1 %>% names() %>% str_match("(sensitivity_wtd)(\\[\\d+\\])") %>% {paste0(.[,2],"_Roche", .[,3])}
	names(dat_Se_Roche_US_age_2) <- dat_Se_Roche_US_age_2 %>% names() %>% str_match("(sensitivity_wtd)(\\[\\d+\\])") %>% {paste0(.[,2],"_Roche", .[,3])}
	names(dat_Se_Roche_US_age_3) <- dat_Se_Roche_US_age_3 %>% names() %>% str_match("(sensitivity_wtd)(\\[\\d+\\])") %>% {paste0(.[,2],"_Roche", .[,3])}
	names(dat_Se_Roche_US_age_4) <- dat_Se_Roche_US_age_4 %>% names() %>% str_match("(sensitivity_wtd)(\\[\\d+\\])") %>% {paste0(.[,2],"_Roche", .[,3])}
	
	## Collate the data from 3 assays
	dat_US_kinetics_age_1 <- cbind(as.matrix(dat_Se_Abbott_US_age_1), as.matrix(dat_Se_VitrosIgG_US_age_1), as.matrix(dat_Se_Roche_US_age_1))
	dat_US_kinetics_age_2 <- cbind(as.matrix(dat_Se_Abbott_US_age_2), as.matrix(dat_Se_VitrosIgG_US_age_2), as.matrix(dat_Se_Roche_US_age_2))
	dat_US_kinetics_age_3 <- cbind(as.matrix(dat_Se_Abbott_US_age_3), as.matrix(dat_Se_VitrosIgG_US_age_3), as.matrix(dat_Se_Roche_US_age_3))
	dat_US_kinetics_age_4 <- cbind(as.matrix(dat_Se_Abbott_US_age_4), as.matrix(dat_Se_VitrosIgG_US_age_4), as.matrix(dat_Se_Roche_US_age_4))
	
	## Estimate the parameters of Se (mean and sd), by assay
	weighted_Se_with_ts_age_1 <- gqs(model_weighted_Se_with_ts, draws=dat_US_kinetics_age_1, data=list(N=N, prop_infections=PROP_INFECTIONS_age_1_0_to_17))
	weighted_Se_with_ts_age_2 <- gqs(model_weighted_Se_with_ts, draws=dat_US_kinetics_age_2, data=list(N=N, prop_infections=PROP_INFECTIONS_age_2_18_to_49))
	weighted_Se_with_ts_age_3 <- gqs(model_weighted_Se_with_ts, draws=dat_US_kinetics_age_3, data=list(N=N, prop_infections=PROP_INFECTIONS_age_3_50_to_64))
	weighted_Se_with_ts_age_4 <- gqs(model_weighted_Se_with_ts, draws=dat_US_kinetics_age_4, data=list(N=N, prop_infections=PROP_INFECTIONS_age_4_65_plus))
	
	## Get raw prevalence
	prev_raw <- dat_sero_US_track$seroprevalence_raw[i]
	
	## Repeat each twice (for F and M)
	MU_SE_Abbott <- c(
		summary(weighted_Se_with_ts_age_1)$summary["sensitivity_wtd_Abbott_with_ts","mean"],
		summary(weighted_Se_with_ts_age_1)$summary["sensitivity_wtd_Abbott_with_ts","mean"],
		summary(weighted_Se_with_ts_age_2)$summary["sensitivity_wtd_Abbott_with_ts","mean"],
		summary(weighted_Se_with_ts_age_2)$summary["sensitivity_wtd_Abbott_with_ts","mean"],
		summary(weighted_Se_with_ts_age_3)$summary["sensitivity_wtd_Abbott_with_ts","mean"],
		summary(weighted_Se_with_ts_age_3)$summary["sensitivity_wtd_Abbott_with_ts","mean"],
		summary(weighted_Se_with_ts_age_4)$summary["sensitivity_wtd_Abbott_with_ts","mean"],
		summary(weighted_Se_with_ts_age_4)$summary["sensitivity_wtd_Abbott_with_ts","mean"])
	
	SIGMA_SE_Abbott <- c(
		summary(weighted_Se_with_ts_age_1)$summary["sensitivity_wtd_Abbott_with_ts","sd"],
		summary(weighted_Se_with_ts_age_1)$summary["sensitivity_wtd_Abbott_with_ts","sd"],
		summary(weighted_Se_with_ts_age_2)$summary["sensitivity_wtd_Abbott_with_ts","sd"],
		summary(weighted_Se_with_ts_age_2)$summary["sensitivity_wtd_Abbott_with_ts","sd"],
		summary(weighted_Se_with_ts_age_3)$summary["sensitivity_wtd_Abbott_with_ts","sd"],
		summary(weighted_Se_with_ts_age_3)$summary["sensitivity_wtd_Abbott_with_ts","sd"],
		summary(weighted_Se_with_ts_age_4)$summary["sensitivity_wtd_Abbott_with_ts","sd"],
		summary(weighted_Se_with_ts_age_4)$summary["sensitivity_wtd_Abbott_with_ts","sd"])
	
	MU_SE_VitrosIgG <- c(
		summary(weighted_Se_with_ts_age_1)$summary["sensitivity_wtd_VitrosIgG_with_ts","mean"],
		summary(weighted_Se_with_ts_age_1)$summary["sensitivity_wtd_VitrosIgG_with_ts","mean"],
		summary(weighted_Se_with_ts_age_2)$summary["sensitivity_wtd_VitrosIgG_with_ts","mean"],
		summary(weighted_Se_with_ts_age_2)$summary["sensitivity_wtd_VitrosIgG_with_ts","mean"],
		summary(weighted_Se_with_ts_age_3)$summary["sensitivity_wtd_VitrosIgG_with_ts","mean"],
		summary(weighted_Se_with_ts_age_3)$summary["sensitivity_wtd_VitrosIgG_with_ts","mean"],
		summary(weighted_Se_with_ts_age_4)$summary["sensitivity_wtd_VitrosIgG_with_ts","mean"],
		summary(weighted_Se_with_ts_age_4)$summary["sensitivity_wtd_VitrosIgG_with_ts","mean"])
	
	SIGMA_SE_VitrosIgG <- c(
		summary(weighted_Se_with_ts_age_1)$summary["sensitivity_wtd_VitrosIgG_with_ts","sd"],
		summary(weighted_Se_with_ts_age_1)$summary["sensitivity_wtd_VitrosIgG_with_ts","sd"],
		summary(weighted_Se_with_ts_age_2)$summary["sensitivity_wtd_VitrosIgG_with_ts","sd"],
		summary(weighted_Se_with_ts_age_2)$summary["sensitivity_wtd_VitrosIgG_with_ts","sd"],
		summary(weighted_Se_with_ts_age_3)$summary["sensitivity_wtd_VitrosIgG_with_ts","sd"],
		summary(weighted_Se_with_ts_age_3)$summary["sensitivity_wtd_VitrosIgG_with_ts","sd"],
		summary(weighted_Se_with_ts_age_4)$summary["sensitivity_wtd_VitrosIgG_with_ts","sd"],
		summary(weighted_Se_with_ts_age_4)$summary["sensitivity_wtd_VitrosIgG_with_ts","sd"])
	
	MU_SE_Roche <- c(
		summary(weighted_Se_with_ts_age_1)$summary["sensitivity_wtd_Roche_with_ts","mean"],
		summary(weighted_Se_with_ts_age_1)$summary["sensitivity_wtd_Roche_with_ts","mean"],
		summary(weighted_Se_with_ts_age_2)$summary["sensitivity_wtd_Roche_with_ts","mean"],
		summary(weighted_Se_with_ts_age_2)$summary["sensitivity_wtd_Roche_with_ts","mean"],
		summary(weighted_Se_with_ts_age_3)$summary["sensitivity_wtd_Roche_with_ts","mean"],
		summary(weighted_Se_with_ts_age_3)$summary["sensitivity_wtd_Roche_with_ts","mean"],
		summary(weighted_Se_with_ts_age_4)$summary["sensitivity_wtd_Roche_with_ts","mean"],
		summary(weighted_Se_with_ts_age_4)$summary["sensitivity_wtd_Roche_with_ts","mean"])
	
	SIGMA_SE_Roche <- c(
		summary(weighted_Se_with_ts_age_1)$summary["sensitivity_wtd_Roche_with_ts","sd"],
		summary(weighted_Se_with_ts_age_1)$summary["sensitivity_wtd_Roche_with_ts","sd"],
		summary(weighted_Se_with_ts_age_2)$summary["sensitivity_wtd_Roche_with_ts","sd"],
		summary(weighted_Se_with_ts_age_2)$summary["sensitivity_wtd_Roche_with_ts","sd"],
		summary(weighted_Se_with_ts_age_3)$summary["sensitivity_wtd_Roche_with_ts","sd"],
		summary(weighted_Se_with_ts_age_3)$summary["sensitivity_wtd_Roche_with_ts","sd"],
		summary(weighted_Se_with_ts_age_4)$summary["sensitivity_wtd_Roche_with_ts","sd"],
		summary(weighted_Se_with_ts_age_4)$summary["sensitivity_wtd_Roche_with_ts","sd"])
	
	## Get the weights by age and sex
	tmp_Abbott %>%
		mutate(prop_age_gender = pop_of_census_division / sum(pop_of_census_division)) %>%
		select(prop_age_gender) %>%
		pull() -> WEIGHTS
	
	if(prev_raw > 0) {
		
		## Fit the model
		## Here: we use the binomial model because we know the numerators and denominators
		fit_Stan <- stan(model_code=calc_seroprev_binomial_3_assays_adjusting_over_n_groups_per_assay,
		data=list(
			n_groups = 8,
			weight = WEIGHTS,
			y_sample_Abbott = NUMERATOR_BY_ASSAY_Abbott,
			n_sample_Abbott = DENOMINATOR_BY_ASSAY_Abbott,
			y_sample_VitrosIgG = NUMERATOR_BY_ASSAY_Ortho,
			n_sample_VitrosIgG = DENOMINATOR_BY_ASSAY_Ortho,
			y_sample_Roche = NUMERATOR_BY_ASSAY_Roche,
			n_sample_Roche = DENOMINATOR_BY_ASSAY_Roche,
			mu_Se_Abbott = MU_SE_Abbott,
			sigma_Se_Abbott = SIGMA_SE_Abbott,
			mu_Se_VitrosIgG = MU_SE_VitrosIgG,
			sigma_Se_VitrosIgG = SIGMA_SE_VitrosIgG,
			mu_Se_Roche = MU_SE_Roche,
			sigma_Se_Roche = SIGMA_SE_Roche,
			# prop_test = PROP_TEST_BY_ASSAY,
			Sp_Abbott = Sp_Abbott,
			Sp_VitrosIgG = Sp_VitrosIgG,
			Sp_Roche = Sp_Roche),
		refresh=0)
		
		## Store
		dat_sero_US_track$adj_sero_median[i] <- summary(fit_Stan)$summary["prev_adj_overall","50%"]
		dat_sero_US_track$adj_sero_lb[i] <- summary(fit_Stan)$summary["prev_adj_overall","2.5%"]
		dat_sero_US_track$adj_sero_ub[i] <- summary(fit_Stan)$summary["prev_adj_overall","97.5%"]
		
	}
	
	else {
		dat_sero_US_track$adj_sero_median[i] <- 0
		dat_sero_US_track$adj_sero_lb[i] <- 0
		dat_sero_US_track$adj_sero_ub[i] <- 0
	}
	
	print(i)
	
}

## Plot
windows()
dat_sero_US_track %>%
	ggplot() +
		geom_pointrange(aes(x=Testing_Period, y=adj_sero_median, ymin = adj_sero_lb, ymax = adj_sero_ub), colour="red") +
		geom_point(aes(x=Testing_Period, y=seroprevalence_raw), colour="black") +
		theme_bw() +
		xlab("") +
		ylab("") +
		facet_grid( ~ Census_Division, scales="fixed") +
		theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Get the slopes of interest
data.frame(intercept=0, slope=1) %>%
	add_row(intercept=0, slope=1.5) %>%
	add_row(intercept=0, slope=2) %>%
	add_row(intercept=0, slope=3) %>%
	mutate(slope_factor = factor(slope, levels=c("3","2","1.5","1"), labels=c("Adj = 3*raw","Adj = 2*raw", "Adj = 1.5*raw","Adj = raw"))) -> dat_slopes

## Scatterplot
windows()
dat_sero_US_track %>%
	ggplot(aes(x=seroprevalence_raw, y=adj_sero_median, colour=Testing_Period)) +
		# geom_abline(dat=dat_slopes, aes(intercept=intercept, slope=slope, linetype=slope_factor), size=0.5) +
		geom_abline(dat=dat_slopes, aes(intercept=intercept, slope=slope, linetype=slope_factor), size=0.5) +
		scale_linetype_manual(name = "Bias", values=rev(c(1,6,5,3))) +
		geom_point() +
		geom_pointrange(aes(ymin=adj_sero_lb, ymax=adj_sero_ub)) +
		theme_bw() +
		scale_x_continuous("Reported seroprevalence", breaks=scales::pretty_breaks(5)) +
		scale_y_continuous("Adjusted seroprevalence", breaks=scales::pretty_breaks(5)) +
		geom_abline(intercept=0, slope=1, colour="black", linetype="dashed") +
		facet_wrap(. ~ Census_Division, scales="fixed")

## ESTIMATE MEAN BY CENSUS DIVISON FOR A SINGLE SURVEY [4]

## CALCULATE WEIGHTED SEROPREV BY AGE
dat_sero_US_track %>%
	filter(Round==4) -> survey_4_weighted_by_age_and_sex

## Re-link them back to states
dat_US_state_counts %>%
	filter(Round==4) %>%
	select(state, census_division) %>%
	left_join(survey_4_weighted_by_age_and_sex, by=c("census_division"="Census_Division")) -> survey_4_weighted_by_age_and_sex_final

survey_4_weighted_by_age_and_sex_final$state <- usdata::abbr2state(survey_4_weighted_by_age_and_sex_final$state)
survey_4_weighted_by_age_and_sex_final$state[49] <- "Puerto Rico"
survey_4_weighted_by_age_and_sex_final$state <- tolower(survey_4_weighted_by_age_and_sex_final$state)

## MAKE A MAP OF THE US
## h/t: https://stackoverflow.com/questions/13757771/relocating-alaska-and-hawaii-on-thematic-map-of-the-usa-with-ggplot2
## https://stackoverflow.com/questions/50499363/how-to-draw-the-outline-of-multiple-us-states-in-r

library(sf)
library(fiftystater)

## What's the largest value to plot
survey_4_weighted_by_age_and_sex_final %>%
	select(-n_positive_raw, -n_tested_raw, -Round) %>%
	select_if(is.numeric) %>%
	max() + 0.01 -> MAX_VALUE_US

sf::st_as_sf(fifty_states, coords = c("long", "lat")) %>% 
  group_by(id, piece) %>% 
  summarize(do_union = FALSE) %>%
  st_cast("POLYGON") %>% 
  ungroup() -> sf_fifty

## Borders (census divisions)
sf_fifty %>%
	left_join(survey_4_weighted_by_age_and_sex_final, by=c("id"="state")) %>%
	filter(id %in% c("hawaii", "alaska") == FALSE) %>%
	group_by(census_division) %>%
	summarise() -> census_division

CD_1 <- as(st_geometry(census_division$geometry[1]), "Spatial")
CD_2 <- as(st_geometry(census_division$geometry[2]), "Spatial")
CD_3 <- as(st_geometry(census_division$geometry[3]), "Spatial")
CD_4 <- as(st_geometry(census_division$geometry[4]), "Spatial")
CD_5 <- as(st_geometry(census_division$geometry[5]), "Spatial")
CD_6 <- as(st_geometry(census_division$geometry[6]), "Spatial")
CD_7 <- as(st_geometry(census_division$geometry[7]), "Spatial")
CD_8 <- as(st_geometry(census_division$geometry[8]), "Spatial")
CD_9 <- as(st_geometry(census_division$geometry[9]), "Spatial")

## Add Puerto Rico
USA <- sf::st_read("Data/tl_2020_us_county.shp")
PR_counties <- USA %>% filter(STATEFP==72)
PR <- as(st_geometry(PR_counties), "Spatial")
PR_outline <- rgeos::gUnaryUnion(PR)
PR_spdf <- broom::tidy(PR_outline)
PR_spdf$adj_sero_median <- survey_4_weighted_by_age_and_sex_final %>% filter(census_division=="South Atlantic & Puerto Rico") %>% slice(1) %>% select(adj_sero_median) %>% pull()

## PLOT
ggplot(as.data.frame(survey_4_weighted_by_age_and_sex_final), aes(map_id=state)) +
  geom_map(aes(fill=adj_sero_median), map = fifty_states, col="black", size=0.5) +
  coord_map() +
  geom_polygon(data=CD_1, aes(x=long, y=lat, group=group, map_id=1), fill=NA, colour="black", size=1, alpha=0) +  
  geom_polygon(data=CD_2, aes(x=long, y=lat, group=group, map_id=1), fill=NA, colour="black", size=1, alpha=0) +  
  geom_polygon(data=CD_3, aes(x=long, y=lat, group=group, map_id=1), fill=NA, colour="black", size=1, alpha=0) +  
  geom_polygon(data=CD_4, aes(x=long, y=lat, group=group, map_id=1), fill=NA, colour="black", size=1, alpha=0) +  
  geom_polygon(data=CD_5, aes(x=long, y=lat, group=group, map_id=1), fill=NA, colour="black", size=1, alpha=0) +  
  geom_polygon(data=CD_6, aes(x=long, y=lat, group=group, map_id=1), fill=NA, colour="black", size=1, alpha=0) +  
  geom_polygon(data=CD_7, aes(x=long, y=lat, group=group, map_id=1), fill=NA, colour="black", size=1, alpha=0) +  
  geom_polygon(data=CD_8, aes(x=long, y=lat, group=group, map_id=1), fill=NA, colour="black", size=1, alpha=0) +  
  geom_polygon(data=CD_9, aes(x=long, y=lat, group=group, map_id=1), fill=NA, colour="black", size=1, alpha=0) +  
  geom_polygon(data=PR_spdf, aes(x=long-9+8, y=lat+7, group=group, map_id=1, fill=adj_sero_median), colour="black", size=0.5) +  
  fifty_states_inset_boxes() +
  geom_rect(aes(xmin=-78+8, xmax=-73+8, ymin=23, ymax=27), colour="black", fill=NA) +
  expand_limits(x = fifty_states$long, y = fifty_states$lat) +
	ggtitle("United States: adjusted seroprevalence (Round 4) - with AS as separate category") +
	theme_bw() +
	theme( ## Box only
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	 axis.ticks = element_blank(),
	axis.text=element_blank(),
	axis.line=element_blank()) +
	xlab("") + ylab("") +
	scale_fill_gradientn("Median", limits=c(0,MAX_VALUE_US), breaks=pretty_breaks(7), colours=brewer.pal(9, "YlOrRd")) -> p_US_adj_4

dat_sero_US_track %>%
	# filter(Round==4) %>%
	rename(`Census Division` = Census_Division) %>%
	mutate(Round = as.factor(Round)) %>%
	ggplot(aes(x=seroprevalence_raw, y=adj_sero_median, colour=`Census Division`, shape=Round)) +
		geom_abline(dat=dat_slopes, aes(intercept=intercept, slope=slope, linetype=slope_factor), size=0.5) +
		# scale_colour_discrete(guide = FALSE) +
		scale_linetype_manual(name = "Bias", values=rev(c(1,6,5,3))) +
		geom_point(size=0.4) +
		geom_pointrange(aes(ymin=adj_sero_lb, ymax=adj_sero_ub)) +
		theme_bw() +
		scale_x_continuous("Raw seroprevalence", breaks=scales::pretty_breaks(5)) +
		scale_y_continuous("Adjusted seroprevalence", breaks=scales::pretty_breaks(5)) -> p_US_scatter_4

p_US_adj_4 / p_US_scatter_4

## Prep to save
dat_sero_US_track_3groups <- dat_sero_US_track
dat_US_state_counts_3groups <- dat_US_state_counts
survey_4_weighted_by_age_and_sex_final_3groups <- survey_4_weighted_by_age_and_sex_final
p_US_adj_4_3groups <- p_US_adj_4
p_US_scatter_4_3groups <- p_US_scatter_4
MAX_VALUE_US_3groups <- MAX_VALUE_US
CD_1_3groups <- CD_1
CD_2_3groups <- CD_2
CD_3_3groups <- CD_3
CD_4_3groups <- CD_4
CD_5_3groups <- CD_5
CD_6_3groups <- CD_6
CD_7_3groups <- CD_7
CD_8_3groups <- CD_8
CD_9_3groups <- CD_9
PR_spdf_3groups <- PR_spdf

## Save for future plotting etc
save(
	dat_sero_US_track_3groups,
	dat_US_state_counts_3groups,
	survey_4_weighted_by_age_and_sex_final_3groups,
	p_US_adj_4_3groups,
	p_US_scatter_4_3groups,
	MAX_VALUE_US_3groups,
	CD_1_3groups, CD_2_3groups, CD_3_3groups,
	CD_4_3groups, CD_5_3groups, CD_6_3groups,
	CD_7_3groups, CD_8_3groups, CD_9_3groups,
	PR_spdf_3groups, dat_slopes,
	file="Results/Results_US_3groups.RData")
