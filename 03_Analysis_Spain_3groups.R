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
FILENAME_Abbott <- "Results/Stan_fit_Abbott_3groups.RData"
load(FILENAME_Abbott)
fit_Stan_Abbott <- fit_Stan

## Read in model to generate new quantities
model_weighted_Se <- stan_model(model_code=calc_wtd_Se_by_time_adjusting_for_severity)

## Fix parameters
prop_hosp <- 0.03821921
prop_AS <- 0.5723913

## Generate severity-weighted sensitivities
weighted_Se_Abbott <- gqs(model_weighted_Se, data=list(prop=c(prop_AS, 1-prop_AS-prop_hosp, prop_hosp)), draws=as.matrix(fit_Stan_Abbott))

## Get the estimates of specificity
Sp_Abbott <- 0.9963

##############################################################################################
## Calculate severity-weighted & epi curve-weighted sensitivity for each assay, and adjust: ##
## 																							                                            ##
## Spain																					                                          ##
##############################################################################################

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

## Remove 3 weeks
SURVEY_DATE_Spain_1 <- ymd("2020-05-04")-weeks(3)
SURVEY_DATE_Spain_2 <- ymd("2020-05-25")-weeks(3)

dat_ts_Spain %>%
	filter(fecha <= SURVEY_DATE_Spain_1) %>%
	mutate(days_before_survey_1 = as.double(SURVEY_DATE_Spain_1 - fecha)) %>%
	group_by(provincia_iso) %>%
	mutate(prop_cases = num_casos / sum(num_casos)) %>%
	ungroup() -> dat_ts_Spain_1_final

dat_ts_Spain %>%
	filter(fecha <= SURVEY_DATE_Spain_2) %>%
	mutate(days_before_survey_2 = as.double(SURVEY_DATE_Spain_2 - fecha)) %>%
	group_by(provincia_iso) %>%
	mutate(prop_cases = num_casos / sum(num_casos)) %>%
	ungroup() -> dat_ts_Spain_2_final

## Bring in serosurvey results
dat_sero_Spain <- read_excel("Data/Seroprevalence_data_Spain.xlsx", sheet=1)
dat_sero_Spain %>%
	rename(Seroprevalence_1 = `Round 1 positive`) %>%
	rename(Seroprevalence_2 = `Round 2 positive`) %>%
	mutate(Seroprevalence_1 = Seroprevalence_1 / 100) %>%
	mutate(Seroprevalence_2 = Seroprevalence_2 / 100) -> dat_sero_Spain

dat_sero_Spain$adj_sero_1_median <- NA
dat_sero_Spain$adj_sero_1_lb <- NA
dat_sero_Spain$adj_sero_1_ub <- NA

dat_sero_Spain$adj_sero_2_median <- NA
dat_sero_Spain$adj_sero_2_lb <- NA
dat_sero_Spain$adj_sero_2_ub <- NA

## Read in model to generate new quantities
## Here: we use the Rogan-Gladen estimate because we don't know the numerators and denominators
model_weighted_Se_with_ts <- stan_model(model_code=calc_wtd_Se_by_time_adjusting_for_severity_and_ts_with_seroprev_RG)

PROVINCES <- unique(dat_ts_Spain_1_final$provincia_iso)

for(i in 1:nrow(dat_sero_Spain)) {
	
	## For each survey
	for(j in 1:2) {
	
	## Get the inputs: time series
	if(j==1) {
		N <- dat_ts_Spain_1_final %>% filter(provincia_iso==PROVINCES[i]) %>% nrow()
		prop_infections <- dat_ts_Spain_1_final %>% filter(provincia_iso==PROVINCES[i]) %>% select(prop_cases) %>% pull()
	}
	if(j==2) {
		N <- dat_ts_Spain_2_final %>% filter(provincia_iso==PROVINCES[i]) %>% nrow()
		prop_infections <- dat_ts_Spain_2_final %>% filter(provincia_iso==PROVINCES[i]) %>% select(prop_cases) %>% pull()
	}
	
	## Get the inputs: sensitivity
	dat_Se_Abbott_Spain <- as.data.frame(weighted_Se_Abbott) %>% select(contains("wtd")) %>% select(1:N) %>% as.matrix()
	
	## Get raw prevalence
	if(j==1) prev_raw <- dat_sero_Spain$Seroprevalence_1[i]
	if(j==2) prev_raw <- dat_sero_Spain$Seroprevalence_2[i]
	
	## Generate new quantities
	weighted_Se_with_ts <- gqs(model_weighted_Se_with_ts, draws=dat_Se_Abbott_Spain, data=list(N=N, prop_infections=prop_infections, prev_raw=prev_raw, Sp=Sp_Abbott))
	
	# windows()
	# print(mcmc_hist(weighted_Se_with_ts) + ggtitle(paste0(dat_sero_Spain$Regioni[i])))
	
	## Store
	if(j==1) {
	dat_sero_Spain$adj_sero_1_median[i] <- summary(weighted_Se_with_ts)$summary["prev_adj","50%"]
	dat_sero_Spain$adj_sero_1_lb[i] <- summary(weighted_Se_with_ts)$summary["prev_adj","2.5%"]
	dat_sero_Spain$adj_sero_1_ub[i] <- summary(weighted_Se_with_ts)$summary["prev_adj","97.5%"]
	}
	if(j==2) {
	dat_sero_Spain$adj_sero_2_median[i] <- summary(weighted_Se_with_ts)$summary["prev_adj","50%"]
	dat_sero_Spain$adj_sero_2_lb[i] <- summary(weighted_Se_with_ts)$summary["prev_adj","2.5%"]
	dat_sero_Spain$adj_sero_2_ub[i] <- summary(weighted_Se_with_ts)$summary["prev_adj","97.5%"]
	}
	}
	
	print(i)
	
}

## Match up names
dat_sero_Spain$Provincia[which(dat_sero_Spain$Provincia=="Bizkaia")] <- "Vizcaya"
dat_sero_Spain$Provincia[which(dat_sero_Spain$Provincia=="Gipuzkoa")] <- "Guipúzcoa"
dat_sero_Spain$Provincia[which(dat_sero_Spain$Provincia=="Tenerife")] <- "Santa Cruz de Tenerife"

## Change negatives to 0
dat_sero_Spain %>%
	mutate(adj_sero_1_median = ifelse(adj_sero_1_median < 0, 0, adj_sero_1_median)) %>%
	mutate(adj_sero_2_median = ifelse(adj_sero_1_median < 0, 0, adj_sero_2_median)) %>%
	mutate(adj_sero_1_lb = ifelse(adj_sero_1_lb < 0, 0, adj_sero_1_lb)) %>%
	mutate(adj_sero_2_lb = ifelse(adj_sero_2_lb < 0, 0, adj_sero_2_lb)) %>%
	mutate(adj_sero_1_ub = ifelse(adj_sero_1_ub < 0, 0, adj_sero_1_ub)) %>%
	mutate(adj_sero_2_ub = ifelse(adj_sero_2_ub < 0, 0, adj_sero_2_ub)) -> dat_sero_Spain

## Plot
Spain <- gadm_sp_loadCountries("ESP", level=2, basefile="./")

Spain$spdf@data <- left_join(Spain$spdf@data, dat_sero_Spain, by=c("NAME_2"="Provincia")) 

Spain_df <- fortify(Spain$spdf, region="NAME_2")
Spain_df <- left_join(Spain_df, Spain$spdf@data, by=c("id"="NAME_2"))

Spain_centroid_names <- aggregate(cbind(long, lat) ~ id, data=Spain_df, FUN=mean)

summary(dat_sero_Spain)

## What's the largest value to plot
dat_sero_Spain %>%
	select_if(is.numeric) %>%
	max() + 0.01 -> MAX_VALUE_SPAIN

ggplot() +
	geom_polygon(data=Spain_df %>% filter(lat >= 35), aes(x=long, y=lat, group=group, fill=Seroprevalence_1), color='black') +
	ggtitle("Raw seroprevalence (Survey 1)") +
	theme_minimal() +
	# geom_text(data=Spain_centroid_names, aes(x=long, y=lat, label=id), size=2.5) +
	scale_fill_gradientn("proportion", limits=c(0,MAX_VALUE_SPAIN), breaks=pretty_breaks(7), colours=brewer.pal(9, "YlOrRd")) -> p_Spain_1_raw

ggplot() +
	geom_polygon(data=Spain_df %>% filter(lat >= 35), aes(x=long, y=lat, group=group, fill=Seroprevalence_2), color='black') +
	ggtitle("Raw seroprevalence (Survey 2)") +
	theme_minimal() +
	# geom_text(data=Spain_centroid_names, aes(x=long, y=lat, label=id), size=2.5) +
	scale_fill_gradientn("proportion", limits=c(0,MAX_VALUE_SPAIN), breaks=pretty_breaks(7), colours=brewer.pal(9, "YlOrRd")) -> p_Spain_2_raw

ggplot() +
	geom_polygon(data=Spain_df %>% filter(lat >= 35), aes(x=long, y=lat, group=group, fill=adj_sero_1_median), color='black') +
	ggtitle("Adjusted seroprevalence 1 (median)") +
	theme_minimal() +
	scale_fill_gradientn("proportion", limits=c(0,MAX_VALUE_SPAIN), breaks=pretty_breaks(7), colours=brewer.pal(9, "YlOrRd")) -> p_Spain_1_adj_median

ggplot() +
	geom_polygon(data=Spain_df %>% filter(lat >= 35), aes(x=long, y=lat, group=group, fill=adj_sero_2_median), color='black') +
	ggtitle("Adjusted seroprevalence 2 (median)") +
	theme_minimal() +
	scale_fill_gradientn("proportion", limits=c(0,MAX_VALUE_SPAIN), breaks=pretty_breaks(7), colours=brewer.pal(9, "YlOrRd")) -> p_Spain_2_adj_median

ggplot() +
	geom_polygon(data=Spain_df %>% filter(lat >= 35), aes(x=long, y=lat, group=group, fill=adj_sero_1_lb), color='black') +
	ggtitle("Adjusted seroprevalence 1 (lb)") +
	theme_minimal() +
	scale_fill_gradientn("proportion", limits=c(0,MAX_VALUE_SPAIN), breaks=pretty_breaks(7), colours=brewer.pal(9, "YlOrRd")) -> p_Spain_1_adj_lb

ggplot() +
	geom_polygon(data=Spain_df %>% filter(lat >= 35), aes(x=long, y=lat, group=group, fill=adj_sero_2_lb), color='black') +
	ggtitle("Adjusted seroprevalence 2 (lb)") +
	theme_minimal() +
	scale_fill_gradientn("proportion", limits=c(0,MAX_VALUE_SPAIN), breaks=pretty_breaks(7), colours=brewer.pal(9, "YlOrRd")) -> p_Spain_2_adj_lb

ggplot() +
	geom_polygon(data=Spain_df %>% filter(lat >= 35), aes(x=long, y=lat, group=group, fill=adj_sero_1_ub), color='black') +
	ggtitle("Adjusted seroprevalence 1 (ub)") +
	theme_minimal() +
	scale_fill_gradientn("proportion", limits=c(0,MAX_VALUE_SPAIN), breaks=pretty_breaks(7), colours=brewer.pal(9, "YlOrRd")) -> p_Spain_1_adj_ub

ggplot() +
	geom_polygon(data=Spain_df %>% filter(lat >= 35), aes(x=long, y=lat, group=group, fill=adj_sero_2_ub), color='black') +
	ggtitle("Adjusted seroprevalence 2 (ub)") +
	theme_minimal() +
	scale_fill_gradientn("proportion", limits=c(0,MAX_VALUE_SPAIN), breaks=pretty_breaks(7), colours=brewer.pal(9, "YlOrRd")) -> p_Spain_2_adj_ub

windows()
(p_Spain_1_raw + p_Spain_1_adj_median) / (p_Spain_1_adj_lb + p_Spain_1_adj_ub) + plot_layout(guides="collect") + plot_annotation(tag_levels='A') + theme(plot.tag=element_text(face='bold'))

windows()
(p_Spain_2_raw + p_Spain_2_adj_median) / (p_Spain_2_adj_lb + p_Spain_2_adj_ub) + plot_layout(guides="collect") + plot_annotation(tag_levels='A') + theme(plot.tag=element_text(face='bold'))

## Get the slopes of interest
data.frame(intercept=0, slope=1) %>%
	add_row(intercept=0, slope=1.5) %>%
	add_row(intercept=0, slope=2) %>%
	add_row(intercept=0, slope=3) %>%
	mutate(slope_factor = factor(slope, levels=c("3","2","1.5","1"), labels=c("Adj = 3*raw", "Adj = 2*raw", "Adj = 1.5*raw","Adj = raw"))) -> dat_slopes

## Scatterplots
dat_sero_Spain %>%
	ggplot(aes(x=Seroprevalence_1, y=adj_sero_1_median)) +
		geom_abline(dat=dat_slopes, aes(intercept=intercept, slope=slope, linetype=slope_factor), size=0.5) +
		geom_point() +
		geom_pointrange(aes(ymin=adj_sero_1_lb, ymax=adj_sero_1_ub)) +
		theme_bw() +
		scale_x_continuous("Reported seroprevalence (survey 1)", breaks=scales::pretty_breaks(10)) +
		scale_y_continuous("Adjusted seroprevalence", breaks=scales::pretty_breaks(10)) -> plot_scatter_1

dat_sero_Spain %>%
	ggplot(aes(x=Seroprevalence_2, y=adj_sero_2_median)) +
		geom_abline(dat=dat_slopes, aes(intercept=intercept, slope=slope, linetype=slope_factor), size=0.5) +
		geom_point() +
		geom_pointrange(aes(ymin=adj_sero_2_lb, ymax=adj_sero_2_ub)) +
		theme_bw() +
		scale_x_continuous("Reported seroprevalence (survey 2)", breaks=scales::pretty_breaks(10)) +
		scale_y_continuous("Adjusted seroprevalence", breaks=scales::pretty_breaks(10)) -> plot_scatter_2

windows()
plot_scatter_1 + plot_scatter_2

## Plot by region and by time
dat_sero_Spain %>%
	group_by(Provincia) %>%
	pivot_longer(Seroprevalence_1:adj_sero_2_ub) %>%
	ungroup() %>%
	filter(str_detect(name, "Seroprevalence")) %>%
	mutate(name = factor(name, levels=c("Seroprevalence_1", "Seroprevalence_2"), labels=c("Round 1", "Round 2"))) -> dat_sero_Spain_reported

dat_sero_Spain %>%
	group_by(Provincia) %>%
	pivot_longer(Seroprevalence_1:adj_sero_2_ub) %>%
	ungroup() %>%
	filter(!str_detect(name, "Seroprevalence")) %>%
	separate(name, c("tmp1","tmp2","name","metric")) %>%
	select(-tmp1, -tmp2) %>%
	group_by(Provincia, name) %>%
	pivot_wider(names_from=metric, values_from=value) %>%
	ungroup() %>%
	mutate(name = factor(name, levels=c("1", "2"), labels=c("Round 1", "Round 2"))) -> dat_sero_Spain_adjusted

#############################
## Putting it all together ##
#############################

## Rename
Spain_df_3groups <- Spain_df
dat_sero_Spain_reported_3groups <- dat_sero_Spain_reported
dat_sero_Spain_adjusted_3groups <- dat_sero_Spain_adjusted
MAX_VALUE_SPAIN_3groups <- MAX_VALUE_SPAIN

ggplot() +
	geom_polygon(data=Spain_df_3groups %>% filter(lat >= 35), aes(x=long, y=lat, group=group, fill=adj_sero_1_median), color='black') +
	ggtitle("Spain: adjusted seroprevalence (Round 1) - with AS as separate category") +
	coord_map() +
	theme_bw() +
	theme( ## Box only
	legend.position = "none",
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	 axis.ticks = element_blank(),
	axis.text=element_blank(),
	axis.line=element_blank()) +
	xlab("") + ylab("") +
	scale_fill_gradientn("Median", limits=c(0,MAX_VALUE_SPAIN_3groups), breaks=pretty_breaks(7), colours=brewer.pal(9, "YlOrRd")) -> p_Spain_1_adj_median_3groups

ggplot() +
	geom_polygon(data=Spain_df_3groups %>% filter(lat >= 35), aes(x=long, y=lat, group=group, fill=adj_sero_2_median), color='black') +
	ggtitle("Spain: adjusted seroprevalence (Round 2) - with AS as separate category") +
	coord_map() +
	theme_bw() +
	theme( ## Box only
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	 axis.ticks = element_blank(),
	axis.text=element_blank(),
	axis.line=element_blank()) +
	xlab("") + ylab("") +
	scale_fill_gradientn("Median", limits=c(0,MAX_VALUE_SPAIN_3groups), breaks=pretty_breaks(7), colours=brewer.pal(9, "YlOrRd")) -> p_Spain_2_adj_median_3groups

dat_sero_Spain_reported_3groups %>%
	rename(reported=value) %>%
	left_join(dat_sero_Spain_adjusted_3groups, by=c("Area","Comunidades Autónomas","Provincia","ISO","name")) %>%
	rename(Round = name) %>%
	ggplot(aes(x=reported, y=median, colour=Round)) +
		geom_abline(dat=dat_slopes, aes(intercept=intercept, slope=slope, linetype=slope_factor), size=0.5) +
		scale_linetype_manual(name = "Bias", values=rev(c(1,6,5,3))) +
		geom_point(size=0.4) +
		geom_pointrange(aes(ymin=lb, ymax=ub)) +
		theme_bw() +
		scale_x_continuous("Raw seroprevalence", breaks=scales::pretty_breaks(5)) +
		scale_y_continuous("Adjusted seroprevalence", breaks=scales::pretty_breaks(5)) -> p_Spain_scatter_3groups

(p_Spain_1_adj_median_3groups / p_Spain_2_adj_median_3groups) | p_Spain_scatter_3groups

## Bring in demographic data
dat_demog_Spain <- read_excel("Data/Demography_Spain_province.xlsx")
dat_demog_Spain %>%
	mutate(Population_pct = Population / sum(Population)) -> dat_demog_Spain

## Estimate single value for country
dat_sero_Spain_adjusted_3groups %>%
	left_join(dat_demog_Spain, by=c("Provincia"="Province")) %>%
	group_by(name) %>%
	summarise(
		median = weighted.mean(median, Population_pct),
		lb = weighted.mean(lb, Population_pct),
		ub = weighted.mean(ub, Population_pct))

## Clean up
dat_sero_Spain_reported_3groups %>%
	rename(reported=value) %>%
	left_join(dat_sero_Spain_adjusted_3groups, by=c("Area","Comunidades Autónomas","Provincia","ISO","name")) %>%
	rename(Round = name) -> dat_sero_Spain_3groups

## Save for future plotting etc
save(Spain_df_3groups, dat_sero_Spain_3groups, dat_sero_Spain_reported_3groups, dat_sero_Spain_adjusted_3groups, MAX_VALUE_SPAIN_3groups, dat_demog_Spain, dat_slopes, p_Spain_1_adj_median_3groups, p_Spain_2_adj_median_3groups, p_Spain_scatter_3groups, file="Results/Results_Spain_3groups.RData")
