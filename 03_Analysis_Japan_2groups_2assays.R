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

## Load the bivariate results
FILENAME <- "Results/Stan_fit_bivariate_2groups.RData"
load(FILENAME)

## Look at traceplots
windows()
mcmc_trace(fit_Stan, pars=vars(contains("lambda"), contains("log_beta_severity"), contains("Sigma"), contains("tau"), contains("theta"))) + theme_bw()

## Look at sensitivity
windows()
as.data.frame(summary(fit_Stan)$summary) %>%
	tibble::rownames_to_column() %>%
	filter(str_detect(rowname, paste0("sensitivity_over_time\\["))) %>%
	mutate(severity_final_class = as.numeric(str_match_all(rowname, "(?<=\\[).+?(?=\\,)"))) %>%
	mutate(day = as.numeric(str_match_all(rowname, "(?<=\\,).+?(?=\\])"))) %>%
	mutate(severity_final_class = factor(severity_final_class, levels=c(1:8), labels=c("Non-hosp (+,+)", "Non-hosp (+,-)", "Non-hosp (-,+)", "Non-hosp (-,-)",
	"Hosp (+,+)", "Hosp (+,-)", "Hosp (-,+)", "Hosp (-,-)"))) %>%
	ggplot() +
		theme_bw() +
		theme(legend.position="none") +
		geom_ribbon(aes(x=day, ymin=`2.5%`, ymax=`97.5%`, group=severity_final_class, fill=severity_final_class), alpha=0.2) +
		geom_line(aes(x=day, y=`50%`, group=severity_final_class, colour=severity_final_class), size=1) +
		scale_x_continuous(breaks=scales::pretty_breaks(4)) +
		scale_y_continuous(breaks=scales::pretty_breaks(5)) +
		xlab("Days - 21") +
		ylab("Sensitivity") +
		ggtitle("(Abbott, Roche)") +
		facet_wrap(~severity_final_class, scales="fixed", ncol=4)

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
dat_sero_Japan_bivariate <- read_excel("Data/Seroprevalence_data_Japan.xlsx", sheet=1)

## Get raw prevalence
dat_sero_Japan_bivariate$bivariate_raw <- dat_sero_Japan_bivariate$Positive_Both / dat_sero_Japan_bivariate$N_tested

dat_sero_Japan_bivariate$bivariate_median <- NA
dat_sero_Japan_bivariate$bivariate_lb <- NA
dat_sero_Japan_bivariate$bivariate_ub <- NA

## Read in model to generate new quantities (bivariate)
calc_wtd_Se_by_time_adjusting_for_severity_2groups_bivariate <-'
data {
	
	real prop[2];
	
}
parameters {
	
	matrix<lower=0.0, upper=1.0>[8,366] sensitivity_over_time;
	
}
generated quantities {
	
	real sensitivity_wtd[4,366];
	
	for(i in 1:366) {
		
		sensitivity_wtd[1,i] = (prop[1]*sensitivity_over_time[1,i]) + (prop[2]*sensitivity_over_time[5,i]);
		sensitivity_wtd[2,i] = (prop[1]*sensitivity_over_time[2,i]) + (prop[2]*sensitivity_over_time[6,i]);
		sensitivity_wtd[3,i] = (prop[1]*sensitivity_over_time[3,i]) + (prop[2]*sensitivity_over_time[7,i]);
		sensitivity_wtd[4,i] = (prop[1]*sensitivity_over_time[4,i]) + (prop[2]*sensitivity_over_time[8,i]);
		
	}
}
'

model_weighted_Se <- stan_model(model_code=calc_wtd_Se_by_time_adjusting_for_severity_2groups_bivariate)

## Generate severity-weighted sensitivities
## This is ordered as (Abbott, Roche): (+,+), (+,-), (-,+), (-,-)
weighted_Se <- gqs(model_weighted_Se, data=list(prop=c(1-prop_hosp, prop_hosp)), draws=as.matrix(fit_Stan))

## Get the inputs: time series
N <- dat_ts_Japan_final %>% nrow()
N_times_4 <- N*4
prop_infections <- dat_ts_Japan_final %>% select(prop_cases) %>% pull()

## Get the inputs: sensitivity
## This is ordered as (Abbott, Roche): (+,+), (+,-), (-,+), (-,-)
dat_Se_Japan <- as.data.frame(weighted_Se) %>% select(contains("wtd")) %>% select(1:N_times_4) %>% as.matrix()

## Read in model to generate new quantities, for the bivariate situation
calc_wtd_Se_by_time_adjusting_for_severity_and_ts_bivariate <-'
data {
	
	int N;
	real prop_infections[N];
	
}
parameters {
	
	real sensitivity_wtd[4,N];
	
}
generated quantities {
	
	real sensitivity_wtd_with_ts[4];
	
	sensitivity_wtd_with_ts[1] = 0.0;
	sensitivity_wtd_with_ts[2] = 0.0;
	sensitivity_wtd_with_ts[3] = 0.0;
	sensitivity_wtd_with_ts[4] = 0.0;
	
	for(i in 1:N) {
		
		sensitivity_wtd_with_ts[1] += sensitivity_wtd[1,i] * prop_infections[N+1-i];
		sensitivity_wtd_with_ts[2] += sensitivity_wtd[2,i] * prop_infections[N+1-i];
		sensitivity_wtd_with_ts[3] += sensitivity_wtd[3,i] * prop_infections[N+1-i];
		sensitivity_wtd_with_ts[4] += sensitivity_wtd[4,i] * prop_infections[N+1-i];
		
	}
	
}
'

model_weighted_Se_with_ts <- stan_model(model_code=calc_wtd_Se_by_time_adjusting_for_severity_and_ts_bivariate)

## Estimate the parameters of Se (mean and sd), by assay
weighted_Se <- gqs(model_weighted_Se_with_ts, draws=dat_Se_Japan, data=list(N=N, prop_infections=prop_infections))

## Get the mu and sigma for each Se
mu_Se_1 <- summary(weighted_Se)$summary["sensitivity_wtd_with_ts[1]","mean"]
mu_Se_2 <- summary(weighted_Se)$summary["sensitivity_wtd_with_ts[2]","mean"]
mu_Se_3 <- summary(weighted_Se)$summary["sensitivity_wtd_with_ts[3]","mean"]
mu_Se_4 <- summary(weighted_Se)$summary["sensitivity_wtd_with_ts[4]","mean"]

sigma_Se_1 <- summary(weighted_Se)$summary["sensitivity_wtd_with_ts[1]","sd"]
sigma_Se_2 <- summary(weighted_Se)$summary["sensitivity_wtd_with_ts[2]","sd"]
sigma_Se_3 <- summary(weighted_Se)$summary["sensitivity_wtd_with_ts[3]","sd"]
sigma_Se_4 <- summary(weighted_Se)$summary["sensitivity_wtd_with_ts[4]","sd"]

## Model adapted from: https://pubmed.ncbi.nlm.nih.gov/10802336/
Stan_model_cond_covariance <- '
data {
  
  // This is ordered as (Abbott, Roche): (+,+), (+,-), (-,+), (-,-)
  real mu_Se[4];
  real sigma_Se[4];
  //real Se[4];
  
  // This is ordered as (Abbott, Roche): (+,+), (+,-), (-,+), (-,-)
  int y_sample_Tokyo[4];
  int y_sample_Osaka[4];
  int y_sample_Miyagi[4];
  int y_sample_Aichi[4];
  int y_sample_Fukuoka[4];
  
  real Sp_Abbott;
  real Sp_Roche;
  
}

parameters {
  
  real<lower=0.0, upper=1.0> prev_Tokyo;
  real<lower=0.0, upper=1.0> prev_Osaka;
  real<lower=0.0, upper=1.0> prev_Miyagi;
  real<lower=0.0, upper=1.0> prev_Aichi;
  real<lower=0.0, upper=1.0> prev_Fukuoka;

  simplex[4] Se;
  
}

transformed parameters {
  
  simplex[4] p_sample_Tokyo;
  simplex[4] p_sample_Osaka;
  simplex[4] p_sample_Miyagi;
  simplex[4] p_sample_Aichi;
  simplex[4] p_sample_Fukuoka;
  
  p_sample_Tokyo[1] = prev_Tokyo*(Se[1]) + (1.0-prev_Tokyo)*(1.0-Sp_Abbott)*(1.0-Sp_Roche);
  p_sample_Tokyo[2] = prev_Tokyo*(Se[2]) + (1.0-prev_Tokyo)*(1.0-Sp_Abbott)*Sp_Roche;
  p_sample_Tokyo[3] = prev_Tokyo*(Se[3]) + (1.0-prev_Tokyo)*Sp_Abbott*(1.0-Sp_Roche);
  p_sample_Tokyo[4] = prev_Tokyo*(Se[4]) + (1.0-prev_Tokyo)*Sp_Abbott*Sp_Roche;
  
  p_sample_Osaka[1] = prev_Osaka*(Se[1]) + (1.0-prev_Osaka)*(1.0-Sp_Abbott)*(1.0-Sp_Roche);
  p_sample_Osaka[2] = prev_Osaka*(Se[2]) + (1.0-prev_Osaka)*(1.0-Sp_Abbott)*Sp_Roche;
  p_sample_Osaka[3] = prev_Osaka*(Se[3]) + (1.0-prev_Osaka)*Sp_Abbott*(1.0-Sp_Roche);
  p_sample_Osaka[4] = prev_Osaka*(Se[4]) + (1.0-prev_Osaka)*Sp_Abbott*Sp_Roche;
  
  p_sample_Miyagi[1] = prev_Miyagi*(Se[1]) + (1.0-prev_Miyagi)*(1.0-Sp_Abbott)*(1.0-Sp_Roche);
  p_sample_Miyagi[2] = prev_Miyagi*(Se[2]) + (1.0-prev_Miyagi)*(1.0-Sp_Abbott)*Sp_Roche;
  p_sample_Miyagi[3] = prev_Miyagi*(Se[3]) + (1.0-prev_Miyagi)*Sp_Abbott*(1.0-Sp_Roche);
  p_sample_Miyagi[4] = prev_Miyagi*(Se[4]) + (1.0-prev_Miyagi)*Sp_Abbott*Sp_Roche;
  
  p_sample_Aichi[1] = prev_Aichi*(Se[1]) + (1.0-prev_Aichi)*(1.0-Sp_Abbott)*(1.0-Sp_Roche);
  p_sample_Aichi[2] = prev_Aichi*(Se[2]) + (1.0-prev_Aichi)*(1.0-Sp_Abbott)*Sp_Roche;
  p_sample_Aichi[3] = prev_Aichi*(Se[3]) + (1.0-prev_Aichi)*Sp_Abbott*(1.0-Sp_Roche);
  p_sample_Aichi[4] = prev_Aichi*(Se[4]) + (1.0-prev_Aichi)*Sp_Abbott*Sp_Roche;
  
  p_sample_Fukuoka[1] = prev_Fukuoka*(Se[1]) + (1.0-prev_Fukuoka)*(1.0-Sp_Abbott)*(1.0-Sp_Roche);
  p_sample_Fukuoka[2] = prev_Fukuoka*(Se[2]) + (1.0-prev_Fukuoka)*(1.0-Sp_Abbott)*Sp_Roche;
  p_sample_Fukuoka[3] = prev_Fukuoka*(Se[3]) + (1.0-prev_Fukuoka)*Sp_Abbott*(1.0-Sp_Roche);
  p_sample_Fukuoka[4] = prev_Fukuoka*(Se[4]) + (1.0-prev_Fukuoka)*Sp_Abbott*Sp_Roche;
}

model {
  
  // Sero-prevalence
  target += multinomial_lpmf(y_sample_Tokyo | p_sample_Tokyo);
  target += multinomial_lpmf(y_sample_Osaka | p_sample_Osaka);
  target += multinomial_lpmf(y_sample_Miyagi | p_sample_Miyagi);
  target += multinomial_lpmf(y_sample_Aichi | p_sample_Aichi);
  target += multinomial_lpmf(y_sample_Fukuoka | p_sample_Fukuoka);

  for(i in 1:4) target += normal_lpdf(Se[i] | mu_Se[i], sigma_Se[i]);
  
}
'

## Fit the model
fit_Stan_cond_covariance <- stan(
	model_code=Stan_model_cond_covariance,
	data=list(
	## This is ordered as (Abbott, Roche): (+,+), (+,-), (-,+), (-,-)
	mu_Se = c(mu_Se_1, mu_Se_2, mu_Se_3, mu_Se_4),
	sigma_Se = c(sigma_Se_1, sigma_Se_2, sigma_Se_3, sigma_Se_4),
	## This is ordered as (Abbott, Roche): (+,+), (+,-), (-,+), (-,-)
	y_sample_Tokyo=c(31,6,29,3333),
	y_sample_Osaka=c(16,5,9,2716),
	y_sample_Miyagi=c(4,8,5,2843),
	y_sample_Aichi=c(16,9,11,2924),
	y_sample_Fukuoka=c(6,6,13,3053),
	## From package inserts
	Sp_Abbott=Sp_Abbott,
	Sp_Roche=Sp_Roche),
	chains=4, iter=5000, thin=5, init=0, control=list(adapt_delta=0.9))

## Store
dat_sero_Japan_bivariate$bivariate_median[1] <- summary(fit_Stan_cond_covariance)$summary["prev_Tokyo","50%"]
dat_sero_Japan_bivariate$bivariate_lb[1] <- summary(fit_Stan_cond_covariance)$summary["prev_Tokyo","2.5%"]
dat_sero_Japan_bivariate$bivariate_ub[1] <- summary(fit_Stan_cond_covariance)$summary["prev_Tokyo","97.5%"]

dat_sero_Japan_bivariate$bivariate_median[2] <- summary(fit_Stan_cond_covariance)$summary["prev_Osaka","50%"]
dat_sero_Japan_bivariate$bivariate_lb[2] <- summary(fit_Stan_cond_covariance)$summary["prev_Osaka","2.5%"]
dat_sero_Japan_bivariate$bivariate_ub[2] <- summary(fit_Stan_cond_covariance)$summary["prev_Osaka","97.5%"]

dat_sero_Japan_bivariate$bivariate_median[3] <- summary(fit_Stan_cond_covariance)$summary["prev_Miyagi","50%"]
dat_sero_Japan_bivariate$bivariate_lb[3] <- summary(fit_Stan_cond_covariance)$summary["prev_Miyagi","2.5%"]
dat_sero_Japan_bivariate$bivariate_ub[3] <- summary(fit_Stan_cond_covariance)$summary["prev_Miyagi","97.5%"]

dat_sero_Japan_bivariate$bivariate_median[4] <- summary(fit_Stan_cond_covariance)$summary["prev_Aichi","50%"]
dat_sero_Japan_bivariate$bivariate_lb[4] <- summary(fit_Stan_cond_covariance)$summary["prev_Aichi","2.5%"]
dat_sero_Japan_bivariate$bivariate_ub[4] <- summary(fit_Stan_cond_covariance)$summary["prev_Aichi","97.5%"]

dat_sero_Japan_bivariate$bivariate_median[5] <- summary(fit_Stan_cond_covariance)$summary["prev_Fukuoka","50%"]
dat_sero_Japan_bivariate$bivariate_lb[5] <- summary(fit_Stan_cond_covariance)$summary["prev_Fukuoka","2.5%"]
dat_sero_Japan_bivariate$bivariate_ub[5] <- summary(fit_Stan_cond_covariance)$summary["prev_Fukuoka","97.5%"]

## Plot
windows()
dat_sero_Japan_bivariate %>%
	pivot_longer(cols=bivariate_raw:bivariate_ub) %>%
	separate(name, sep="_", into=c("Assay","metric")) %>%
	filter(metric=="raw") %>%
	mutate(metric = "raw") %>%
	mutate(Assay = "Bivariate") %>%
	ggplot() +
		geom_pointrange(data=dat_sero_Japan_bivariate %>% mutate(Assay="Bivariate") %>% mutate(metric="adjusted"), aes(x=Assay, y=bivariate_median, ymin=bivariate_lb, ymax=bivariate_ub, colour=metric)) +
		# geom_point(data=dat_sero_Japan_bivariate %>% mutate(Assay="bivariate") %>% mutate(metric="raw (Abbott +)"), aes(x=Assay, y=Positive_Abbott/N_tested, colour=metric), size=2.5) +
		# geom_point(data=dat_sero_Japan_bivariate %>% mutate(Assay="bivariate") %>% mutate(metric="raw (Roche +)"), aes(x=Assay, y=Positive_Roche/N_tested, colour=metric), size=2.5) +
		geom_point(aes(x=Assay, y=value, colour=metric), size=2.5) +
		scale_colour_manual(values=c("raw (Abbott +)"="black", "raw (Roche +)"="black", "raw"="black", "adjusted"="red"), breaks=c("adjusted","raw (Roche +)","raw (Abbott +)","raw")) +
		theme_bw() +
		xlab("") +
		ylab("Seroprevalence") +
		facet_wrap(.~Prefecture, scales="fixed", ncol=5)

## Prep to save
dat_sero_Japan_bivariate_2groups <- dat_sero_Japan_bivariate

# Save for future plotting etc
save(
	dat_sero_Japan_bivariate_2groups,
	file="Results/Results_Japan_2groups_2assays.RData")
