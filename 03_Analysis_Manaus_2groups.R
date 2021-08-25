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
FILENAME_Abbott <- "Results/Stan_fit_Abbott_2groups.RData"
load(FILENAME_Abbott)
fit_Stan_Abbott <- fit_Stan

## Read in model to generate new quantities
## All this is done below
model_weighted_Se <- stan_model(model_code=calc_wtd_Se_by_time_adjusting_for_severity_2groups)

## Get the estimates of specificity
Sp_Abbott <- 0.9963

## Bring in serosurvey results (continuous, monthly)
dat_sero_Manaus_raw <- read_delim("Data/Seroprevalence_data_Manaus.csv", delim=";")

## Fix parameters, by age
## Age groups in Manaus:
## < 15 15-24 25-34 35-44 45-54 55-64 65-70
## Correspondingly:
## Take the mean between 2 consecutive points

severity_by_age <- read_excel("Data/Clinical_fraction_and_severity_by_age.xlsx")

ggplot(severity_by_age, aes(x=Age_group, y=Clinical_fraction_mean)) + geom_point()

## Get the mean between 2 consecutive points
sev_age_Manaus <- data.frame(Age_group=names(table(dat_sero_Manaus_raw$age_range)))
sev_age_Manaus$Clinical_fraction_mean <- NA
sev_age_Manaus$Prob_hosp_given_inf_Salje <- NA

for(i in 1:nrow(sev_age_Manaus)) {
	
	sev_age_Manaus$Clinical_fraction_mean[i] <- mean(c(severity_by_age$Clinical_fraction_mean[i], severity_by_age$Clinical_fraction_mean[i+1]))
	sev_age_Manaus$Prob_hosp_given_inf_Salje[i] <- mean(c(severity_by_age$Prob_hosp_given_inf_Salje[i], severity_by_age$Prob_hosp_given_inf_Salje[i+1]))
	
}

bind_rows(severity_by_age %>% select(Age_group, Clinical_fraction_mean, Prob_hosp_given_inf_Salje) %>% mutate(group="original"), sev_age_Manaus %>% mutate(group="Manaus")) %>%
	mutate(Age_group = factor(Age_group, levels=c("0_to_9","< 15","10_to_19","15-24","20_to_29","25-34","30_to_39","35-44","40_to_49","45-54","50_to_59","55-64","60_to_69","65-70","70_to_79","80_and_up"))) -> dat_to_plot_severity_by_age

dat_to_plot_severity_by_age %>%
	ggplot(aes(x=Age_group, y=Clinical_fraction_mean, colour=group)) +
		geom_point(size=4) +
		theme_bw() -> p1

dat_to_plot_severity_by_age %>%
	ggplot(aes(x=Age_group, y=Prob_hosp_given_inf_Salje, colour=group)) +
		geom_point(size=4) +
		theme_bw() -> p2

p1 / p2

## Generate severity-weighted sensitivities, by age

sev_age_Manaus$prop_AS <- 1-sev_age_Manaus$Clinical_fraction_mean
sev_age_Manaus$prop_hosp <- sev_age_Manaus$Prob_hosp_given_inf_Salje

weighted_Se_Abbott_1_under_15 <- gqs(model_weighted_Se, data=list(prop=c(1-sev_age_Manaus$prop_hosp[1], sev_age_Manaus$prop_hosp[1])), draws=as.matrix(fit_Stan_Abbott))

weighted_Se_Abbott_2_15_to_24 <- gqs(model_weighted_Se, data=list(prop=c(1-sev_age_Manaus$prop_hosp[2], sev_age_Manaus$prop_hosp[2])), draws=as.matrix(fit_Stan_Abbott))

weighted_Se_Abbott_3_25_to_34 <- gqs(model_weighted_Se, data=list(prop=c(1-sev_age_Manaus$prop_hosp[3], sev_age_Manaus$prop_hosp[3])), draws=as.matrix(fit_Stan_Abbott))

weighted_Se_Abbott_4_35_to_44 <- gqs(model_weighted_Se, data=list(prop=c(1-sev_age_Manaus$prop_hosp[4], sev_age_Manaus$prop_hosp[4])), draws=as.matrix(fit_Stan_Abbott))

weighted_Se_Abbott_5_45_to_54 <- gqs(model_weighted_Se, data=list(prop=c(1-sev_age_Manaus$prop_hosp[5], sev_age_Manaus$prop_hosp[5])), draws=as.matrix(fit_Stan_Abbott))

weighted_Se_Abbott_6_55_to_64 <- gqs(model_weighted_Se, data=list(prop=c(1-sev_age_Manaus$prop_hosp[6], sev_age_Manaus$prop_hosp[6])), draws=as.matrix(fit_Stan_Abbott))

weighted_Se_Abbott_7_65_to_70 <- gqs(model_weighted_Se, data=list(prop=c(1-sev_age_Manaus$prop_hosp[7], sev_age_Manaus$prop_hosp[7])), draws=as.matrix(fit_Stan_Abbott))

##############################################################################################
## Calculate severity-weighted & epi curve-weighted sensitivity for each assay, and adjust: ##
## 																						                                            	##
## Manaus, Brazil																	                                      		##
##############################################################################################

## Plot
windows()
dat_sero_Manaus_raw %>%
	mutate(clean_date = ymd(paste0(donation_year, "-", donation_month, "-01"))) %>%
	group_by(clean_date, gender, age_range) %>%
	summarise(
		n_tested = n()) %>%
	ggplot(aes(x=clean_date, y=n_tested, group=age_range)) +
		geom_line(aes(colour=age_range)) +
		geom_point(aes(colour=age_range), size=2) +
		facet_wrap(. ~ gender, scales="fixed", ncol=1) +
		theme_bw() +
		xlab("") +
		scale_x_date(date_labels="%b %Y", date_breaks="1 month") +
		theme(axis.text.x=element_text(angle=45, hjust=1)) +
		ggtitle("Sample size")

dat_sero_Manaus_raw %>%
	mutate(clean_date = ymd(paste0(donation_year, "-", donation_month, "-15"))) %>%
	mutate(seropos = ifelse(abbott_signal_to_cutoff >= 1.4, 1, 0)) %>%
	# group_by(clean_date) %>%
	group_by(clean_date, gender, age_range) %>%
	summarise(
		raw_prev = sum(seropos) / n(),
		y_sample = sum(seropos),
		n_sample = n()) %>%
	ungroup() -> dat_sero_Manaus

# dat_sero_Manaus$adj_sero_median <- NA
# dat_sero_Manaus$adj_sero_lb <- NA
# dat_sero_Manaus$adj_sero_ub <- NA

# dat_sero_Manaus$weighted_Se_median <- NA
# dat_sero_Manaus$weighted_Se_lb <- NA
# dat_sero_Manaus$weighted_Se_ub <- NA

## Bring in epi curve (hospitalizations)
read_excel("Data/TimeSeries_data_Manaus_hospitalizations.xlsx", sheet=1) %>%
	mutate(date=as.Date(date)) %>%
	arrange(date) %>%
	filter(date >= ymd("2020-03-01")) %>%
	complete(date=seq.Date(min(date), max(date), by="day"), fill=list(hosp=0)) -> dat_ts_Manaus_hosp_raw

## Bring in epi curve from EpiNow2
load("Results/EpiNow2_fit_Manaus.RData")

backcalc_using_reported_hosp_Bi["summarised"] %>%
	as.data.frame() %>%
	filter(summarised.variable=="infections") %>%
	droplevels() %>%
	mutate(summarised.date = ymd(summarised.date)) %>%
	mutate(param="Bi") -> data_frame_Manaus_EpiNow2_Bi

backcalc_using_reported_hosp_Linton["summarised"] %>%
	as.data.frame() %>%
	filter(summarised.variable=="infections") %>%
	droplevels() %>%
	mutate(summarised.date = ymd(summarised.date)) %>%
	mutate(param="Linton") -> data_frame_Manaus_EpiNow2_Linton

## Combine the data (medians only)
bind_rows(data_frame_Manaus_EpiNow2_Bi, data_frame_Manaus_EpiNow2_Linton) %>%
	rename(date = summarised.date) %>%
	rename(onsets = summarised.median) %>%
	mutate(model = "EpiNow2") %>%
	select(date, onsets, model, param) -> dat_ts_Manaus_2_ways

## Pick a time series
dat_ts_Manaus_2_ways %>%
	filter(model=="EpiNow2" & param=="Linton") %>%
	droplevels() -> dat_ts_Manaus

## Monthly surveys, remove 3 weeks
SURVEY_DATE_Manaus_1 <- ymd("2020-03-15")-weeks(3)
SURVEY_DATE_Manaus_2 <- ymd("2020-04-15")-weeks(3)
SURVEY_DATE_Manaus_3 <- ymd("2020-05-15")-weeks(3)
SURVEY_DATE_Manaus_4 <- ymd("2020-06-15")-weeks(3)
SURVEY_DATE_Manaus_5 <- ymd("2020-07-15")-weeks(3)
SURVEY_DATE_Manaus_6 <- ymd("2020-08-15")-weeks(3)
SURVEY_DATE_Manaus_7 <- ymd("2020-09-15")-weeks(3)
SURVEY_DATE_Manaus_8 <- ymd("2020-10-15")-weeks(3)
SURVEY_DATE_Manaus_9 <- ymd("2020-11-15")-weeks(3)
SURVEY_DATE_Manaus_10 <- ymd("2020-12-15")-weeks(3)
SURVEY_DATE_Manaus_11 <- ymd("2021-01-15")-weeks(3)

dat_ts_Manaus %>%
	filter(date <= SURVEY_DATE_Manaus_1) %>%
	mutate(days_before_survey_1 = as.double(SURVEY_DATE_Manaus_1 - date)) %>%
	mutate(prop_cases = onsets / sum(onsets)) -> dat_ts_Manaus_1_final

dat_ts_Manaus %>%
	filter(date <= SURVEY_DATE_Manaus_2) %>%
	mutate(days_before_survey_2 = as.double(SURVEY_DATE_Manaus_2 - date)) %>%
	mutate(prop_cases = onsets / sum(onsets)) -> dat_ts_Manaus_2_final

dat_ts_Manaus %>%
	filter(date <= SURVEY_DATE_Manaus_3) %>%
	mutate(days_before_survey_3 = as.double(SURVEY_DATE_Manaus_3 - date)) %>%
	mutate(prop_cases = onsets / sum(onsets)) -> dat_ts_Manaus_3_final

dat_ts_Manaus %>%
	filter(date <= SURVEY_DATE_Manaus_4) %>%
	mutate(days_before_survey_4 = as.double(SURVEY_DATE_Manaus_4 - date)) %>%
	mutate(prop_cases = onsets / sum(onsets)) -> dat_ts_Manaus_4_final

dat_ts_Manaus %>%
	filter(date <= SURVEY_DATE_Manaus_5) %>%
	mutate(days_before_survey_5 = as.double(SURVEY_DATE_Manaus_5 - date)) %>%
	mutate(prop_cases = onsets / sum(onsets)) -> dat_ts_Manaus_5_final

dat_ts_Manaus %>%
	filter(date <= SURVEY_DATE_Manaus_6) %>%
	mutate(days_before_survey_6 = as.double(SURVEY_DATE_Manaus_6 - date)) %>%
	mutate(prop_cases = onsets / sum(onsets)) -> dat_ts_Manaus_6_final

dat_ts_Manaus %>%
	filter(date <= SURVEY_DATE_Manaus_7) %>%
	mutate(days_before_survey_7 = as.double(SURVEY_DATE_Manaus_7 - date)) %>%
	mutate(prop_cases = onsets / sum(onsets)) -> dat_ts_Manaus_7_final

dat_ts_Manaus %>%
	filter(date <= SURVEY_DATE_Manaus_8) %>%
	mutate(days_before_survey_8 = as.double(SURVEY_DATE_Manaus_8 - date)) %>%
	mutate(prop_cases = onsets / sum(onsets)) -> dat_ts_Manaus_8_final

dat_ts_Manaus %>%
	filter(date <= SURVEY_DATE_Manaus_9) %>%
	mutate(days_before_survey_9 = as.double(SURVEY_DATE_Manaus_9 - date)) %>%
	mutate(prop_cases = onsets / sum(onsets)) -> dat_ts_Manaus_9_final

dat_ts_Manaus %>%
	filter(date <= SURVEY_DATE_Manaus_10) %>%
	mutate(days_before_survey_10 = as.double(SURVEY_DATE_Manaus_10 - date)) %>%
	mutate(prop_cases = onsets / sum(onsets)) -> dat_ts_Manaus_10_final

dat_ts_Manaus %>%
	filter(date <= SURVEY_DATE_Manaus_11) %>%
	mutate(days_before_survey_11 = as.double(SURVEY_DATE_Manaus_11 - date)) %>%
	mutate(prop_cases = onsets / sum(onsets)) -> dat_ts_Manaus_11_final

## Read in model to generate new quantities
## This is just giving us the mean and SD of the weighted (single) sensitivity
model_weighted_Se_with_ts <- stan_model(model_code=calc_wtd_Se_by_time_adjusting_for_severity_and_ts)

## REMOVING MARCH SURVEY
dat_sero_Manaus %>%
	filter(clean_date != "2020-03-15") -> dat_sero_Manaus

## REMOVE YOUNGEST AGE BIN
dat_sero_Manaus %>%
	filter(age_range != "< 15") -> dat_sero_Manaus

## NEED TO COMPLETE VARS
dat_sero_Manaus %>%
	arrange(clean_date, age_range, gender) %>%
	complete(clean_date, gender, age_range, fill=list(raw_prev=NA, y_sample=0, n_sample=0)) %>%
	arrange(clean_date, age_range, gender) -> dat_sero_Manaus

## Bring in demographic data
dat_pop_by_age <- read_excel("Data/Demography_Manaus_age_sex.xlsx")

dat_pop_by_age$age_range <- c("< 15","< 15","< 15","15-24","15-24","25-34","25-34","35-44","35-44","45-54","45-54","55-64","55-64","65-70")

dat_pop_by_age %>%
	pivot_longer(cols=Male_N:Female_N, names_to="gender", values_to="N") %>%
	mutate(gender = ifelse(gender=="Male_N", "M", "F")) %>%
	filter(age_range != "< 15") %>%
	group_by(age_range, gender) %>%
	summarise(
		N_age_gender = sum(N)) %>%
	ungroup() %>%
	mutate(prop_age_gender = N_age_gender / sum(N_age_gender)) -> dat_pop_by_age_clean

dat_sero_Manaus %>%
	left_join(dat_pop_by_age_clean, by=c("age_range","gender")) -> dat_sero_Manaus

## Aggregate over each time point
dat_sero_Manaus %>%
	group_by(clean_date) %>%
	summarise(raw_prev = sum(y_sample) / sum(n_sample)) -> dat_sero_Manaus_AGG

dat_sero_Manaus_AGG$adj_sero_median <- NA
dat_sero_Manaus_AGG$adj_sero_lb <- NA
dat_sero_Manaus_AGG$adj_sero_ub <- NA

for(i in 1:nrow(dat_sero_Manaus_AGG)) {
	
	## Get the inputs: time series
	if(dat_sero_Manaus_AGG$clean_date[i]=="2020-03-15") dat_FOR_THIS <- dat_ts_Manaus_1_final
	if(dat_sero_Manaus_AGG$clean_date[i]=="2020-04-15") dat_FOR_THIS <- dat_ts_Manaus_2_final
	if(dat_sero_Manaus_AGG$clean_date[i]=="2020-05-15") dat_FOR_THIS <- dat_ts_Manaus_3_final
	if(dat_sero_Manaus_AGG$clean_date[i]=="2020-06-15") dat_FOR_THIS <- dat_ts_Manaus_4_final
	if(dat_sero_Manaus_AGG$clean_date[i]=="2020-07-15") dat_FOR_THIS <- dat_ts_Manaus_5_final
	if(dat_sero_Manaus_AGG$clean_date[i]=="2020-08-15") dat_FOR_THIS <- dat_ts_Manaus_6_final
	if(dat_sero_Manaus_AGG$clean_date[i]=="2020-09-15") dat_FOR_THIS <- dat_ts_Manaus_7_final
	if(dat_sero_Manaus_AGG$clean_date[i]=="2020-10-15") dat_FOR_THIS <- dat_ts_Manaus_8_final
	if(dat_sero_Manaus_AGG$clean_date[i]=="2020-11-15") dat_FOR_THIS <- dat_ts_Manaus_9_final
	if(dat_sero_Manaus_AGG$clean_date[i]=="2020-12-15") dat_FOR_THIS <- dat_ts_Manaus_10_final
	if(dat_sero_Manaus_AGG$clean_date[i]=="2021-01-15") dat_FOR_THIS <- dat_ts_Manaus_11_final
	
	N <- dat_FOR_THIS %>% nrow()
	prop_infections <- dat_FOR_THIS %>% select(prop_cases) %>% pull()
	
	## Get the inputs: sensitivity
	## THIS DEPENDS ON AGE
	# dat_Se_Abbott_Manaus_age_1 <- as.data.frame(weighted_Se_Abbott_1_under_15) %>% select(contains("wtd")) %>% select(1:N) %>% as.matrix()
	dat_Se_Abbott_Manaus_age_2 <- as.data.frame(weighted_Se_Abbott_2_15_to_24) %>% select(contains("wtd")) %>% select(1:N) %>% as.matrix()
	dat_Se_Abbott_Manaus_age_3 <- as.data.frame(weighted_Se_Abbott_3_25_to_34) %>% select(contains("wtd")) %>% select(1:N) %>% as.matrix()
	dat_Se_Abbott_Manaus_age_4 <- as.data.frame(weighted_Se_Abbott_4_35_to_44) %>% select(contains("wtd")) %>% select(1:N) %>% as.matrix()
	dat_Se_Abbott_Manaus_age_5 <- as.data.frame(weighted_Se_Abbott_5_45_to_54) %>% select(contains("wtd")) %>% select(1:N) %>% as.matrix()
	dat_Se_Abbott_Manaus_age_6 <- as.data.frame(weighted_Se_Abbott_6_55_to_64) %>% select(contains("wtd")) %>% select(1:N) %>% as.matrix()
	dat_Se_Abbott_Manaus_age_7 <- as.data.frame(weighted_Se_Abbott_7_65_to_70) %>% select(contains("wtd")) %>% select(1:N) %>% as.matrix()
	
	## Estimate the parameters of Se (mean and sd)
	# weighted_Se_with_ts_age_1 <- gqs(model_weighted_Se_with_ts, draws=dat_Se_Abbott_Manaus_age_1, data=list(N=N, prop_infections=prop_infections))
	weighted_Se_with_ts_age_2 <- gqs(model_weighted_Se_with_ts, draws=dat_Se_Abbott_Manaus_age_2, data=list(N=N, prop_infections=prop_infections))
	weighted_Se_with_ts_age_3 <- gqs(model_weighted_Se_with_ts, draws=dat_Se_Abbott_Manaus_age_3, data=list(N=N, prop_infections=prop_infections))
	weighted_Se_with_ts_age_4 <- gqs(model_weighted_Se_with_ts, draws=dat_Se_Abbott_Manaus_age_4, data=list(N=N, prop_infections=prop_infections))
	weighted_Se_with_ts_age_5 <- gqs(model_weighted_Se_with_ts, draws=dat_Se_Abbott_Manaus_age_5, data=list(N=N, prop_infections=prop_infections))
	weighted_Se_with_ts_age_6 <- gqs(model_weighted_Se_with_ts, draws=dat_Se_Abbott_Manaus_age_6, data=list(N=N, prop_infections=prop_infections))
	weighted_Se_with_ts_age_7 <- gqs(model_weighted_Se_with_ts, draws=dat_Se_Abbott_Manaus_age_7, data=list(N=N, prop_infections=prop_infections))
	
	# windows()
	# print(mcmc_hist(weighted_Se_with_ts) + ggtitle(paste0(dat_sero_Manaus$clean_date[i])))
	
	## Store
	# dat_sero_Manaus$weighted_Se_median[i] <- summary(weighted_Se_with_ts)$summary["sensitivity_wtd_with_ts","50%"]
	# dat_sero_Manaus$weighted_Se_lb[i] <- summary(weighted_Se_with_ts)$summary["sensitivity_wtd_with_ts","2.5%"]
	# dat_sero_Manaus$weighted_Se_ub[i] <- summary(weighted_Se_with_ts)$summary["sensitivity_wtd_with_ts","97.5%"]
	
	## Get raw prevalence
	prev_raw <- dat_sero_Manaus_AGG$raw_prev[i]
	
	MU_SE <- c(
		summary(weighted_Se_with_ts_age_2)$summary["sensitivity_wtd_with_ts","mean"],
		summary(weighted_Se_with_ts_age_2)$summary["sensitivity_wtd_with_ts","mean"],
		summary(weighted_Se_with_ts_age_3)$summary["sensitivity_wtd_with_ts","mean"],
		summary(weighted_Se_with_ts_age_3)$summary["sensitivity_wtd_with_ts","mean"],
		summary(weighted_Se_with_ts_age_4)$summary["sensitivity_wtd_with_ts","mean"],
		summary(weighted_Se_with_ts_age_4)$summary["sensitivity_wtd_with_ts","mean"],
		summary(weighted_Se_with_ts_age_5)$summary["sensitivity_wtd_with_ts","mean"],
		summary(weighted_Se_with_ts_age_5)$summary["sensitivity_wtd_with_ts","mean"],
		summary(weighted_Se_with_ts_age_6)$summary["sensitivity_wtd_with_ts","mean"],
		summary(weighted_Se_with_ts_age_6)$summary["sensitivity_wtd_with_ts","mean"],
		summary(weighted_Se_with_ts_age_7)$summary["sensitivity_wtd_with_ts","mean"],
		summary(weighted_Se_with_ts_age_7)$summary["sensitivity_wtd_with_ts","mean"])
	
	SIGMA_SE <- c(
		summary(weighted_Se_with_ts_age_2)$summary["sensitivity_wtd_with_ts","sd"],
		summary(weighted_Se_with_ts_age_2)$summary["sensitivity_wtd_with_ts","sd"],
		summary(weighted_Se_with_ts_age_3)$summary["sensitivity_wtd_with_ts","sd"],
		summary(weighted_Se_with_ts_age_3)$summary["sensitivity_wtd_with_ts","sd"],
		summary(weighted_Se_with_ts_age_4)$summary["sensitivity_wtd_with_ts","sd"],
		summary(weighted_Se_with_ts_age_4)$summary["sensitivity_wtd_with_ts","sd"],
		summary(weighted_Se_with_ts_age_5)$summary["sensitivity_wtd_with_ts","sd"],
		summary(weighted_Se_with_ts_age_5)$summary["sensitivity_wtd_with_ts","sd"],
		summary(weighted_Se_with_ts_age_6)$summary["sensitivity_wtd_with_ts","sd"],
		summary(weighted_Se_with_ts_age_6)$summary["sensitivity_wtd_with_ts","sd"],
		summary(weighted_Se_with_ts_age_7)$summary["sensitivity_wtd_with_ts","sd"],
		summary(weighted_Se_with_ts_age_7)$summary["sensitivity_wtd_with_ts","sd"])
	
	if(prev_raw > 0) {
		
		## Fit the model
		## Here: we use the binomial model because we know the numerators and denominators
		fit_Stan <- stan(model_code=calc_seroprev_binomial_adjusting_over_n_groups,
		data=list(
			n_groups = as.integer(12),
			weight = dat_sero_Manaus$prop_age_gender[c(((i*12)-11):(i*12))],
			y_sample = dat_sero_Manaus$y_sample[c(((i*12)-11):(i*12))],
			n_sample = dat_sero_Manaus$n_sample[c(((i*12)-11):(i*12))],
			mu_Se = MU_SE,
			sigma_Se = SIGMA_SE,
			Sp = Sp_Abbott),
		refresh=0)
		
		## Store
		dat_sero_Manaus_AGG$adj_sero_median[i] <- summary(fit_Stan)$summary["prev_adj_overall","50%"]
		dat_sero_Manaus_AGG$adj_sero_lb[i] <- summary(fit_Stan)$summary["prev_adj_overall","2.5%"]
		dat_sero_Manaus_AGG$adj_sero_ub[i] <- summary(fit_Stan)$summary["prev_adj_overall","97.5%"]
		
	}
	
	else {
		dat_sero_Manaus_AGG$adj_sero_median[i] <- 0
		dat_sero_Manaus_AGG$adj_sero_lb[i] <- 0
		dat_sero_Manaus_AGG$adj_sero_ub[i] <- 0
	}
	
	print(i)
	
}

## Plot
dat_sero_Manaus_AGG %>%
	ggplot() +
		theme_bw() +
		theme(legend.position="none") +
		geom_pointrange(aes(x=clean_date-days(14), y=adj_sero_median, ymin=adj_sero_lb, ymax=adj_sero_ub), size=0.5) +
		scale_x_date(date_labels="%b %Y", date_breaks="1 month") +
		xlab("") +
		ggtitle("Manaus, Brazil: adjusted seroprevalence over time, across ages 15-70") +
		ylab("Adjusted seroprevalence") -> p_Manaus_ts

as.Date_origin <- function(x){
  format(as.Date(x, origin = '1970-01-01'), '%b %Y')
}

## Get the slopes of interest
data.frame(intercept=0, slope=1) %>%
	add_row(intercept=0, slope=1.5) %>%
	add_row(intercept=0, slope=2) %>%
	add_row(intercept=0, slope=3) %>%
	mutate(slope_factor = factor(slope, levels=c("3","2","1.5","1"), labels=c("Adj = 3*raw","Adj = 2*raw", "Adj = 1.5*raw","Adj = raw"))) -> dat_slopes

dat_sero_Manaus_AGG %>%
	ggplot(aes(x=raw_prev, y=adj_sero_median, colour=as.factor(clean_date - days(14)))) +
		geom_abline(dat=dat_slopes, aes(intercept=intercept, slope=slope, linetype=slope_factor), size=0.5) +
		# scale_colour_discrete(guide = FALSE) +
		scale_linetype_manual(name = "Bias", values=rev(c(1,6,5,3))) +
		geom_point(size=4) +
		geom_linerange(aes(ymin=adj_sero_lb, ymax=adj_sero_ub), size=0.75) +
		theme_bw() +
		scale_colour_brewer("Month", palette="Spectral", labels=as.Date_origin) + 
		scale_x_continuous("Raw seroprevalence", breaks=scales::pretty_breaks(5)) +
		scale_y_continuous("Adjusted seroprevalence", breaks=scales::pretty_breaks(5)) -> p_Manaus_scatter

p_Manaus_ts / p_Manaus_scatter

## Prep to save
dat_sero_Manaus_2groups <- dat_sero_Manaus
dat_sero_Manaus_AGG_2groups <- dat_sero_Manaus_AGG
p_Manaus_ts_2groups <- p_Manaus_ts
p_Manaus_scatter_2groups <- p_Manaus_scatter

## Save for future plotting etc
save(
	dat_sero_Manaus_2groups,
	dat_sero_Manaus_AGG_2groups,
	p_Manaus_ts_2groups,
	p_Manaus_scatter_2groups,
	dat_slopes,
	as.Date_origin,
	file="Results/Results_Manaus_2groups.RData")
