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

## Read in model to generate new quantities
model_weighted_Se <- stan_model(model_code=calc_wtd_Se_by_time_adjusting_for_severity_2groups)

## Fix parameters
prop_hosp <- 0.04389209

## Generate severity-weighted sensitivities
weighted_Se_Abbott <- gqs(model_weighted_Se, data=list(prop=c(1-prop_hosp, prop_hosp)), draws=as.matrix(fit_Stan_Abbott))

## Get the estimates of specificity
Sp_Abbott <- 0.9963

##############################################################################################
## Calculate severity-weighted & epi curve-weighted sensitivity for each assay, and adjust: ##
## 																							                                            ##
## Italy																					                                          ##
##############################################################################################

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

## Remove 3 weeks
SURVEY_DATE_Italy <- ymd("2020-06-19")-weeks(3)

dat_ts_Italy %>%
	arrange(DATA_INIZIO_SINTOMI) %>%
	filter(DATA_INIZIO_SINTOMI <= SURVEY_DATE_Italy) %>%
	mutate(days_before_survey = as.double(SURVEY_DATE_Italy -  DATA_INIZIO_SINTOMI)) %>%
	mutate(prop_cases = CASI / sum(CASI)) -> dat_ts_Italy_final

## Bring in serosurvey results
dat_sero_Italy <- read_excel("Data/Seroprevalence_data_Italy.xlsx", sheet=1)

dat_sero_Italy <- dat_sero_Italy %>% filter(Regioni != "00_ITALIA")

dat_sero_Italy$Seroprevalence <- dat_sero_Italy$Seroprevalence / 100

dat_sero_Italy$adj_sero_median <- NA
dat_sero_Italy$adj_sero_lb <- NA
dat_sero_Italy$adj_sero_ub <- NA

dat_sero_Italy$RATIO_median <- NA
dat_sero_Italy$RATIO_lb <- NA
dat_sero_Italy$RATIO_ub <- NA

## Read in model to generate new quantities
## Here: we use the Rogan-Gladen estimate because we don't know the numerators and denominators
model_weighted_Se_with_ts <- stan_model(model_code=calc_wtd_Se_by_time_adjusting_for_severity_and_ts_with_seroprev_RG)

## Get the inputs: time series
N <- dat_ts_Italy_final %>% nrow()
prop_infections <- dat_ts_Italy_final %>% select(prop_cases) %>% pull()

## Get the inputs: sensitivity
dat_Se_Abbott_Italy <- as.data.frame(weighted_Se_Abbott) %>% select(contains("wtd")) %>% select(1:N) %>% as.matrix()

for(i in 1:nrow(dat_sero_Italy)) {
	
	## Get raw prevalence
	prev_raw <- dat_sero_Italy$Seroprevalence[i]
	
	## Generate new quantities
	weighted_Se_with_ts <- gqs(model_weighted_Se_with_ts, draws=dat_Se_Abbott_Italy, data=list(N=N, prop_infections=prop_infections, prev_raw=prev_raw, Sp=Sp_Abbott))
	
	# windows()
	# print(mcmc_hist(weighted_Se_with_ts) + ggtitle(paste0(dat_sero_Italy$Regioni[i])))
	
	## Store
	dat_sero_Italy$adj_sero_median[i] <- summary(weighted_Se_with_ts)$summary["prev_adj","50%"]
	dat_sero_Italy$adj_sero_lb[i] <- summary(weighted_Se_with_ts)$summary["prev_adj","2.5%"]
	dat_sero_Italy$adj_sero_ub[i] <- summary(weighted_Se_with_ts)$summary["prev_adj","97.5%"]
	
	dat_sero_Italy$RATIO_median[i] <- summary(weighted_Se_with_ts)$summary["prev_RATIO","50%"]
	dat_sero_Italy$RATIO_lb[i] <- summary(weighted_Se_with_ts)$summary["prev_RATIO","2.5%"]
	dat_sero_Italy$RATIO_ub[i] <- summary(weighted_Se_with_ts)$summary["prev_RATIO","97.5%"]
	
}

## Match up names
dat_sero_Italy$Regioni[which(dat_sero_Italy$Regioni=="Puglia")] <- "Apulia"
dat_sero_Italy$Regioni[which(dat_sero_Italy$Regioni=="Sicilia")] <- "Sicily"
dat_sero_Italy$Regioni[which(dat_sero_Italy$Regioni=="Trento")] <- "Trentino-Alto Adige"
dat_sero_Italy$Regioni[which(dat_sero_Italy$Regioni=="Valle d'Aosta/Vallée d'Aoste")] <- "Valle d'Aosta"

## Change negatives to 0
dat_sero_Italy %>%
	mutate(adj_sero_median = ifelse(adj_sero_median < 0, 0, adj_sero_median)) %>%
	mutate(adj_sero_lb = ifelse(adj_sero_lb < 0, 0, adj_sero_lb)) %>%
	mutate(adj_sero_ub = ifelse(adj_sero_ub < 0, 0, adj_sero_ub)) -> dat_sero_Italy

## Plot
Italy <- gadm_sp_loadCountries("ITA", level=1, basefile="./")

Italy$spdf@data <- left_join(Italy$spdf@data, dat_sero_Italy, by=c("NAME_1"="Regioni")) 

Italy_df <- fortify(Italy$spdf, region="NAME_1")
Italy_df <- left_join(Italy_df, Italy$spdf@data, by=c("id"="NAME_1"))

Italy_centroid_names <- aggregate(cbind(long, lat) ~ id, data=Italy_df, FUN=mean)

summary(dat_sero_Italy)

## What's the largest seroprevalence value to plot
dat_sero_Italy %>%
	select(-Population_N, -Population_pct) %>%
	select_if(is.numeric) %>%
	max() + 0.01 -> MAX_VALUE_ITALY

ggplot() +
	geom_polygon(data=Italy_df, aes(x=long, y=lat, group=group, fill=Seroprevalence), color='black') +
	ggtitle("Raw seroprevalence") +
	theme_minimal() +
	# geom_text(data=Italy_centroid_names, aes(x=long, y=lat, label=id), size=2.5) +
	scale_fill_gradientn("proportion", limits=c(0,MAX_VALUE_ITALY), breaks=pretty_breaks(7), colours=brewer.pal(9, "YlOrRd")) -> p_Italy_raw

ggplot() +
	geom_polygon(data=Italy_df, aes(x=long, y=lat, group=group, fill=adj_sero_median), color='black') +
	ggtitle("Adjusted seroprevalence (median)") +
	theme_minimal() +
	scale_fill_gradientn("proportion", limits=c(0,MAX_VALUE_ITALY), breaks=pretty_breaks(7), colours=brewer.pal(9, "YlOrRd")) -> p_Italy_adj_median

ggplot() +
	geom_polygon(data=Italy_df, aes(x=long, y=lat, group=group, fill=adj_sero_lb), color='black') +
	ggtitle("Adjusted seroprevalence (lb)") +
	theme_minimal() +
	scale_fill_gradientn("proportion", limits=c(0,MAX_VALUE_ITALY), breaks=pretty_breaks(7), colours=brewer.pal(9, "YlOrRd")) -> p_Italy_adj_lb

ggplot() +
	geom_polygon(data=Italy_df, aes(x=long, y=lat, group=group, fill=adj_sero_ub), color='black') +
	ggtitle("Adjusted seroprevalence (ub)") +
	theme_minimal() +
	scale_fill_gradientn("proportion", limits=c(0,MAX_VALUE_ITALY), breaks=pretty_breaks(7), colours=brewer.pal(9, "YlOrRd")) -> p_Italy_adj_ub

windows()
(p_Italy_raw + p_Italy_adj_median) / (p_Italy_adj_lb + p_Italy_adj_ub) + plot_layout(guides="collect") + plot_annotation(tag_levels='A') + theme(plot.tag=element_text(face='bold'))

## Get the slopes of interest
data.frame(intercept=0, slope=1) %>%
	add_row(intercept=0, slope=1.5) %>%
	add_row(intercept=0, slope=2) %>%
	add_row(intercept=0, slope=3) %>%
	mutate(slope_factor = factor(slope, levels=c("3","2","1.5","1"), labels=c("Adj = 3*raw","Adj = 2*raw", "Adj = 1.5*raw","Adj = raw"))) -> dat_slopes

## Scatterplot
windows()
dat_sero_Italy %>%
	ggplot(aes(x=Seroprevalence, y=adj_sero_median)) +
		geom_abline(dat=dat_slopes %>% rename(Bias = slope_factor), aes(intercept=intercept, slope=slope, colour=Bias), linetype="dashed", size=1) +
		geom_point(size=3) +
		geom_pointrange(aes(ymin=adj_sero_lb, ymax=adj_sero_ub)) +
		theme_bw() +
		scale_x_continuous("Reported seroprevalence", breaks=scales::pretty_breaks(5)) +
		scale_y_continuous("Adjusted seroprevalence", breaks=scales::pretty_breaks(6)) -> p_Italy_scatter

#############################
## Putting it all together ##
#############################

## Rename
Italy_df_2groups <- Italy_df
dat_sero_Italy_2groups <- dat_sero_Italy
MAX_VALUE_ITALY_2groups <- MAX_VALUE_ITALY

ggplot() +
	geom_polygon(data=Italy_df_2groups, aes(x=long, y=lat, group=group, fill=adj_sero_median), color='black') +
	ggtitle("Italy: adjusted seroprevalence") +
	coord_map() +
	theme_bw() +
	theme( ## Box only
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	 axis.ticks = element_blank(),
	axis.text=element_blank(),
	axis.line=element_blank()) +
	xlab("") + ylab("") +
	scale_fill_gradientn("Median", limits=c(0,MAX_VALUE_ITALY_2groups), breaks=pretty_breaks(7), colours=brewer.pal(9, "YlOrRd")) -> p_Italy_adj_median_2groups

dat_sero_Italy_2groups %>%
	ggplot(aes(x=Seroprevalence, y=adj_sero_median)) +
		geom_abline(dat=dat_slopes, aes(intercept=intercept, slope=slope, linetype=slope_factor), size=0.5) +
		scale_linetype_manual(name = "Bias", values=rev(c(1,6,5,3))) +
		geom_point(size=2) +
		geom_pointrange(aes(ymin=adj_sero_lb, ymax=adj_sero_ub)) +
		theme_bw() +
		scale_x_continuous("Raw seroprevalence", breaks=scales::pretty_breaks(5)) +
		scale_y_continuous("Adjusted seroprevalence", breaks=scales::pretty_breaks(5)) -> p_Italy_scatter_2groups

p_Italy_adj_median_2groups / p_Italy_scatter_2groups

## Estimate single value for country
dat_sero_Italy_2groups %>%
	summarise(
		median = weighted.mean(adj_sero_median, Population_pct),
		lb = weighted.mean(adj_sero_lb, Population_pct),
		ub = weighted.mean(adj_sero_ub, Population_pct))

## Save for future plotting etc
save(Italy_df_2groups, dat_sero_Italy_2groups, MAX_VALUE_ITALY_2groups, dat_slopes, p_Italy_adj_median_2groups, p_Italy_scatter_2groups, file="Results/Results_Italy_2groups.RData")
