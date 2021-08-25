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

###########
## Italy ##
###########

## Bring in epi curve (symptom onsets)
read_excel("Data/TimeSeries_data_Italy_covid_19-iss.xlsx", sheet=4) %>%
	select(-iss_date) %>%
	mutate(
		DATA_INIZIO_SINTOMI = dmy(DATA_INIZIO_SINTOMI),
		CASI = as.numeric(CASI)) %>%
	arrange(DATA_INIZIO_SINTOMI) %>%
	slice(-n()) %>%
	complete(DATA_INIZIO_SINTOMI=seq.Date(min(DATA_INIZIO_SINTOMI), max(DATA_INIZIO_SINTOMI), by="day"), fill=list(CASI=0)) -> dat_ts_Italy

ex <- interval(ymd("2020-05-25"), ymd("2020-07-15"))
int_start(ex) + (int_end(ex) - int_start(ex))/2

## Clean up
SURVEY_DATE_Italy <- ymd("2020-06-19")

dat_ts_Italy %>%
	arrange(DATA_INIZIO_SINTOMI) %>%
	filter(DATA_INIZIO_SINTOMI <= SURVEY_DATE_Italy) %>%
	mutate(days_before_survey = as.double(SURVEY_DATE_Italy -  DATA_INIZIO_SINTOMI)) %>%
	mutate(prop_cases = CASI / sum(CASI)) -> dat_ts_Italy_final

dat_ts_Italy %>%
	arrange(DATA_INIZIO_SINTOMI) %>%
	filter(DATA_INIZIO_SINTOMI <= ymd("2020-07-15")) %>%
	ggplot(aes(x=DATA_INIZIO_SINTOMI, y=CASI)) +
		geom_segment(aes(x=as.Date("2020-05-25"), xend=as.Date("2020-07-15"), y=0, yend=0), colour="cornflowerblue", alpha=0.01, size=1) +
		geom_point(aes(x=ymd("2020-06-19"), y=0), colour="blue", size=2.5, alpha=0.5) +
		geom_line(size=0.2) +
		ggtitle("Italy: reported data") +
		theme_bw() +
		xlab("") +
		ylab("Daily symptom onsets") +
		scale_x_date(date_labels="%b %Y") -> plot.Italy

###########
## Spain ##
###########

## Bring in epi curve (symptom onsets)
read_csv("Data/TimeSeries_data_Spain_datos_provincias.csv") %>%
	arrange(fecha, provincia_iso) %>%
	complete(fecha=seq.Date(min(fecha), max(fecha), by="day"), fill=list(num_casos=0)) %>%
	select(fecha, provincia_iso, num_casos) -> dat_ts_Spain

## Need to rename...
dat_ts_Spain$provincia_iso[which(is.na(dat_ts_Spain$provincia_iso))] <- "NAV"
dat_ts_Spain %>% droplevels() -> dat_ts_Spain

ex <- interval(ymd("2020-04-27"), ymd("2020-05-11"))
int_start(ex) + (int_end(ex) - int_start(ex))/2

ex <- interval(ymd("2020-05-18"), ymd("2020-06-01"))
int_start(ex) + (int_end(ex) - int_start(ex))/2

dat_ts_Spain %>%
	left_join(read_excel("Data/Seroprevalence_data_Spain.xlsx"), by=c("provincia_iso"="ISO")) %>%
	filter(fecha <= ymd("2020-06-01")) %>%
	droplevels() %>%
	group_by(fecha) %>%
	summarise(num_casos_tot = sum(num_casos)) %>%
	ungroup() %>%
	ggplot(aes(x=fecha, y=num_casos_tot)) +
		geom_segment(aes(x=as.Date("2020-04-27"), xend=as.Date("2020-05-11"), y=0, yend=0), colour="cornflowerblue", alpha=0.01, size=1) +
		geom_segment(aes(x=as.Date("2020-05-18"), xend=as.Date("2020-06-01"), y=0, yend=0), colour="cornflowerblue", alpha=0.01, size=1) +
		geom_point(aes(x=ymd("2020-05-04"), y=0), colour="blue", size=2.5, alpha=0.5) +
		geom_point(aes(x=ymd("2020-05-25"), y=0), colour="blue", size=2.5, alpha=0.5) +
		geom_line(size=0.2) +
		ggtitle("Spain: reported data") +
		theme_bw() +
		xlab("") +
		ylab("Daily symptom onsets") +
		scale_x_date(date_labels="%b %Y") -> plot.Spain

########
## US ##
########

## Read in merged EpiNow2 data
load("Results/EpiNow2_fit_US.RData")

## Get reported deaths per census division
REPORTED_DEATHS %>%
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
	rename(Time = date) %>%
	group_by(census_division, Time) %>%
	summarise(
		sum_reported_deaths = sum(confirm)) %>%
	group_by(census_division) %>%
	mutate(
		prop_reported_deaths = sum_reported_deaths / sum(sum_reported_deaths)) -> US_deaths_report_census_division

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
	rename(
		onsets = summarised.median,
		onsets_lb = summarised.lower_95,
		onsets_ub = summarised.upper_95) %>%
	droplevels() %>%
	select(Time, onsets, onsets_lb, onsets_ub, state, census_division) -> data_frame_all_states_NYT_DEATHS

ex <- interval(ymd("2020-07-27"), ymd("2020-08-13"))
int_start(ex) + (int_end(ex) - int_start(ex))/2

ex <- interval(ymd("2020-08-10"), ymd("2020-08-27"))
int_start(ex) + (int_end(ex) - int_start(ex))/2

ex <- interval(ymd("2020-08-24"), ymd("2020-09-10"))
int_start(ex) + (int_end(ex) - int_start(ex))/2

ex <- interval(ymd("2020-09-07"), ymd("2020-09-24"))
int_start(ex) + (int_end(ex) - int_start(ex))/2

data_frame_all_states_NYT_DEATHS %>%
	filter(Time <= ymd("2020-09-24")) %>%
	droplevels() %>%
	group_by(census_division, Time) %>%
	summarise(
		sum_onsets = sum(onsets),
		sum_onsets_lb = sum(onsets_lb),
		sum_onsets_ub = sum(onsets_ub)) %>%
	group_by(census_division) %>%
	mutate(
		prop = sum_onsets / sum(sum_onsets),
		prop_lb = sum_onsets_lb / sum(sum_onsets_lb),
		prop_ub = sum_onsets_ub / sum(sum_onsets_ub)) %>%
	ungroup() %>%
	ggplot(aes(x=Time, y=sum_onsets, group=census_division)) +
		## Reported data
		geom_line(data=US_deaths_report_census_division %>% filter(Time <= ymd("2020-09-24")), aes(x=Time, y=sum_reported_deaths, group=census_division), size=0.2) +
		geom_segment(aes(x=as.Date("2020-07-27"), xend=as.Date("2020-08-13"), y=0, yend=0), colour="cornflowerblue", alpha=0.01, size=1) +
		geom_segment(aes(x=as.Date("2020-08-10"), xend=as.Date("2020-08-27"), y=0, yend=0), colour="cornflowerblue", alpha=0.01, size=1) +
		geom_segment(aes(x=as.Date("2020-08-24"), xend=as.Date("2020-09-10"), y=0, yend=0), colour="cornflowerblue", alpha=0.01, size=1) +
		geom_segment(aes(x=as.Date("2020-09-07"), xend=as.Date("2020-09-24"), y=0, yend=0), colour="cornflowerblue", alpha=0.01, size=1) +
		geom_point(aes(x=ymd("2020-08-04"), y=0), colour="blue", size=2.5, alpha=0.5) +
		geom_point(aes(x=ymd("2020-08-18"), y=0), colour="blue", size=2.5, alpha=0.5) +
		geom_point(aes(x=ymd("2020-09-01"), y=0), colour="blue", size=2.5, alpha=0.5) +
		geom_point(aes(x=ymd("2020-09-15"), y=0), colour="blue", size=2.5, alpha=0.5) +
		geom_ribbon(aes(ymin=sum_onsets_lb, ymax=sum_onsets_ub), fill="firebrick", colour=NA, alpha=0.25) +
		geom_line(colour="firebrick", size=1) +
		ggtitle("United States: reconstructed from death reports") +
		theme_bw() +
		theme(legend.position = "none") +
		xlab("") +
		facet_wrap(.~census_division, ncol=3, scales="free_y") +
		# ylab("Daily symptom onsets (normalized)") +
		ylab("Daily symptom onsets") +
		guides(colour=guide_legend(title=NULL, ncol=3)) +
		scale_x_date(date_labels="%b %Y", date_breaks="2 months") -> plot.US_by_CensusDiv

############
## Manaus ##
############

## Bring in epi curve from EpiNow2
load("Results/EpiNow2_fit_Manaus.RData")

backcalc_using_reported_hosp_Linton["summarised"] %>%
	as.data.frame() %>%
	filter(summarised.variable=="infections") %>%
	droplevels() %>%
	mutate(summarised.date = ymd(summarised.date)) %>%
	mutate(param="Linton") -> data_frame_Manaus_EpiNow2_Linton

data_frame_Manaus_EpiNow2_Linton %>%
	ggplot(aes(x=summarised.date, y=summarised.median, group=param)) +
		## Reported data
		geom_line(data=HOSP_DATA, aes(x=date, y=confirm, group=1), size=0.2) +
		geom_ribbon(aes(ymin=summarised.lower_95, ymax=summarised.upper_95), fill="firebrick", colour=NA, alpha=0.25) +
		geom_line(colour="firebrick", size=1) +
		scale_x_date(date_labels="%b %Y", date_breaks="2 months") +
		geom_segment(aes(x=as.Date("2020-04-01"), xend=as.Date("2021-01-31"), y=0, yend=0), colour="cornflowerblue", alpha=0.01, size=1) +
		# geom_point(aes(x=ymd("2020-03-15"), y=0), colour="blue", size=2.5, alpha=0.5) +
		geom_point(aes(x=ymd("2020-04-15"), y=0), colour="blue", size=2.5, alpha=0.5) +
		geom_point(aes(x=ymd("2020-05-15"), y=0), colour="blue", size=2.5, alpha=0.5) +
		geom_point(aes(x=ymd("2020-06-15"), y=0), colour="blue", size=2.5, alpha=0.5) +
		geom_point(aes(x=ymd("2020-07-15"), y=0), colour="blue", size=2.5, alpha=0.5) +
		geom_point(aes(x=ymd("2020-08-15"), y=0), colour="blue", size=2.5, alpha=0.5) +
		geom_point(aes(x=ymd("2020-09-15"), y=0), colour="blue", size=2.5, alpha=0.5) +
		geom_point(aes(x=ymd("2020-10-15"), y=0), colour="blue", size=2.5, alpha=0.5) +
		geom_point(aes(x=ymd("2020-11-15"), y=0), colour="blue", size=2.5, alpha=0.5) +
		geom_point(aes(x=ymd("2020-12-15"), y=0), colour="blue", size=2.5, alpha=0.5) +
		geom_point(aes(x=ymd("2021-01-15"), y=0), colour="blue", size=2.5, alpha=0.5) +
		xlab("") +
		ylab("Daily symptom onsets") +
		theme_bw() +
		theme(legend.position = "none") +
		ggtitle("Manaus, Brazil: reconstructed from hospitalization reports") -> plot.Manaus

############
## Japan ##
############

## Bring in epi curve from EpiNow2
load("Results/EpiNow2_fit_Japan.RData")

backcalc_using_reported_deaths["summarised"] %>%
	as.data.frame() %>%
	filter(summarised.variable=="infections") %>%
	droplevels() %>%
	mutate(summarised.date = ymd(summarised.date)) %>%
	mutate(param="Linton") -> data_frame_Japan_EpiNow2_Linton

data_frame_Japan_EpiNow2_Linton %>%
filter(summarised.date <= ymd("2021-01-01")) %>%
	ggplot(aes(x=summarised.date, y=summarised.median, group=param)) +
		## Reported data
		geom_line(data=DEATH_DATA %>% filter(date <= ymd("2021-01-01")), aes(x=date, y=confirm, group=1), size=0.2) +
		geom_ribbon(aes(ymin=summarised.lower_95, ymax=summarised.upper_95), fill="firebrick", colour=NA, alpha=0.25) +
		geom_line(colour="firebrick", size=1) +
		scale_x_date(date_labels="%b %Y", date_breaks="2 months") +
		geom_segment(aes(x=as.Date("2020-12-14"), xend=as.Date("2020-12-25"), y=0, yend=0), colour="cornflowerblue", alpha=0.01, size=1) +
		geom_point(aes(x=ymd("2020-12-19"), y=0), colour="blue", size=2.5, alpha=0.5) +
		xlab("") +
		ylab("Daily symptom onsets") +
		theme_bw() +
		theme(legend.position = "none") +
		ggtitle("Japan: reconstructed from death reports") -> plot.Japan

## Put it together
(plot.Italy + plot.Spain) / plot.US_by_CensusDiv / (plot.Manaus + plot.Japan) + plot_layout(heights=c(1,3,1)) + plot_annotation(tag_levels = 'A')
