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
library(loo)

## Shared params for plotting
max_days <- 140
n_breaks <- 6

#######################
## Abbott [2 groups] ##
#######################

FILENAME_Abbott <- "Results/Stan_fit_Abbott_2groups.RData"
load(FILENAME_Abbott)

fit_Stan.Abbott <- fit_Stan
dat_for_Stan_final.Abbott <- dat_for_Stan_final
rm(fit_Stan, dat_for_Stan_final)

## Get the max days by severity
## They're all the same, so just save one
dat_for_Stan_final.Abbott %>%
	group_by(severity_final_class) %>%
	summarise(max_days = max(time_since)) %>%
	ungroup() %>%
	mutate(severity_final = factor(severity_final_class, levels=c(1,2), labels=c("Not hospitalized", "Hospitalized"))) %>%
	mutate(severity_final_class = severity_final) -> tmp_vertical

# dat_for_Stan_final.Roche %>%
	# group_by(severity_final_class) %>%
	# summarise(max_days = max(time_since))

# dat_for_Stan_final.VitrosIgG %>%
	# group_by(severity_final_class) %>%
	# summarise(max_days = max(time_since))

## Plot kinetics
dat_for_Stan_final.Abbott %>%
	group_by(severity_final_class, patient_ID_index, time_since) %>%
	summarise(Mean_log_SC = mean(log(SC_index))) %>%
	ungroup() %>%
	mutate(severity_final = factor(severity_final_class, levels=c(1,2), labels=c("Not hospitalized", "Hospitalized"))) %>%
	ggplot(aes(x=time_since, y=Mean_log_SC, group=patient_ID_index, fill=severity_final)) +
		geom_hline(yintercept=log(1.4), colour="black", linetype="dashed") +
		geom_vline(data=tmp_vertical, aes(xintercept=max_days), colour="purple") +
		geom_line(color="black") +
		geom_point(size=2, shape=21) +
		scale_fill_manual(values=c("cornflowerblue","red")) +
		theme_bw() +
		theme(legend.position="none") +
		facet_wrap(~severity_final) +
		scale_x_continuous(limits=c(0,max_days), breaks=scales::pretty_breaks(n_breaks)) +
		ggtitle("Abbott ARCHITECT") +
		xlab("Days") +
		ylab("log S/C value") -> plot_kinetics.Abbott

## Plot sensitivity
as.data.frame(summary(fit_Stan.Abbott)$summary) %>%
	tibble::rownames_to_column() %>%
	filter(str_detect(rowname, paste0("sensitivity_over_time\\["))) %>%
	mutate(severity_final_class = as.numeric(str_match_all(rowname, "(?<=\\[).+?(?=\\,)"))) %>%
	mutate(day = as.numeric(str_match_all(rowname, "(?<=\\,).+?(?=\\])"))) %>%
	mutate(severity_final_class = factor(severity_final_class, levels=c(1,2), labels=c("Not hospitalized", "Hospitalized"))) %>%
	# filter(day <= 140) %>%
	ggplot() +
		theme_bw() +
		theme(legend.position="none") +
		geom_vline(data=tmp_vertical, aes(xintercept=max_days), colour="purple") +
		geom_ribbon(aes(x=day, ymin=`2.5%`, ymax=`97.5%`, group=severity_final_class, fill=severity_final_class), alpha=0.2) +
		geom_line(aes(x=day, y=`50%`, group=severity_final_class, colour=severity_final_class), size=1) +
		scale_x_continuous(breaks=scales::pretty_breaks(n_breaks)) +
		scale_y_continuous(breaks=scales::pretty_breaks(5), limits=c(0,1)) +
		scale_fill_manual(values=c("Not hospitalized"="cornflowerblue", "Hospitalized"="red")) +
		scale_colour_manual(values=c("Not hospitalized"="cornflowerblue", "Hospitalized"="red")) +
		xlab("Days") +
		ylab("Sensitivity") +
		facet_wrap(~severity_final_class, scales="fixed", ncol=2) -> plot_Sens.Abbott

######################
## Roche [2 groups] ##
######################

FILENAME_Roche <- "Results/Stan_fit_Roche_2groups.RData"
load(FILENAME_Roche)

fit_Stan.Roche <- fit_Stan
dat_for_Stan_final.Roche <- dat_for_Stan_final
rm(fit_Stan, dat_for_Stan_final)

## Plot kinetics
dat_for_Stan_final.Roche %>%
	group_by(severity_final_class, patient_ID, time_since) %>%
	summarise(Mean_log_SC = mean(log(SC_index))) %>%
	ungroup() %>%
	mutate(severity_final = factor(severity_final_class, levels=c(1,2), labels=c("Not hospitalized", "Hospitalized"))) %>%
	ggplot(aes(x=time_since, y=Mean_log_SC, group=patient_ID, fill=severity_final)) +
		geom_hline(yintercept=log(1), colour="black", linetype="dashed") +
		geom_vline(data=tmp_vertical, aes(xintercept=max_days), colour="purple") +
		geom_line(color="black") +
		geom_point(size=2, shape=21) +
		scale_fill_manual(values=c("cornflowerblue","red")) +
		theme_bw() +
		theme(legend.position="none") +
		facet_wrap(~severity_final) +
		scale_x_continuous(limits=c(0,max_days), breaks=scales::pretty_breaks(n_breaks)) +
		ggtitle("Roche Elecsys") +
		xlab("Days") +
		ylab("log S/C value") -> plot_kinetics.Roche

## Plot sensitivity
as.data.frame(summary(fit_Stan.Roche)$summary) %>%
	tibble::rownames_to_column() %>%
	filter(str_detect(rowname, paste0("sensitivity_over_time\\["))) %>%
	mutate(severity_final_class = as.numeric(str_match_all(rowname, "(?<=\\[).+?(?=\\,)"))) %>%
	mutate(day = as.numeric(str_match_all(rowname, "(?<=\\,).+?(?=\\])"))) %>%
	mutate(severity_final_class = factor(severity_final_class, levels=c(1,2), labels=c("Not hospitalized", "Hospitalized"))) %>%
	# filter(day <= 140) %>%
	ggplot() +
		theme_bw() +
		theme(legend.position="none") +
		geom_vline(data=tmp_vertical, aes(xintercept=max_days), colour="purple") +
		geom_ribbon(aes(x=day, ymin=`2.5%`, ymax=`97.5%`, group=severity_final_class, fill=severity_final_class), alpha=0.2) +
		geom_line(aes(x=day, y=`50%`, group=severity_final_class, colour=severity_final_class), size=1) +
		scale_x_continuous(breaks=scales::pretty_breaks(n_breaks)) +
		scale_y_continuous(breaks=scales::pretty_breaks(5), limits=c(0,1)) +
		scale_fill_manual(values=c("Not hospitalized"="cornflowerblue", "Hospitalized"="red")) +
		scale_colour_manual(values=c("Not hospitalized"="cornflowerblue", "Hospitalized"="red")) +
		xlab("Days") +
		ylab("Sensitivity") +
		facet_wrap(~severity_final_class, scales="fixed", ncol=2) -> plot_Sens.Roche

###########################
## Vitros IgG [2 groups] ##
###########################

FILENAME_VitrosIgG <- "Results/Stan_fit_VitrosIgG_2groups.RData"
load(FILENAME_VitrosIgG)

fit_Stan.VitrosIgG <- fit_Stan
dat_for_Stan_final.VitrosIgG <- dat_for_Stan_final
rm(fit_Stan, dat_for_Stan_final)

## Plot kinetics
dat_for_Stan_final.VitrosIgG %>%
	group_by(severity_final_class, patient_ID, time_since) %>%
	summarise(Mean_log_SC = mean(log(SC_index))) %>%
	ungroup() %>%
	mutate(severity_final = factor(severity_final_class, levels=c(1,2), labels=c("Not hospitalized", "Hospitalized"))) %>%
	ggplot(aes(x=time_since, y=Mean_log_SC, group=patient_ID, fill=severity_final)) +
		geom_hline(yintercept=log(1), colour="black", linetype="dashed") +
		geom_vline(data=tmp_vertical, aes(xintercept=max_days), colour="purple") +
		geom_line(color="black") +
		geom_point(size=2, shape=21) +
		scale_fill_manual(values=c("cornflowerblue","red")) +
		theme_bw() +
		theme(legend.position="none") +
		facet_wrap(~severity_final) +
		scale_x_continuous(limits=c(0,max_days), breaks=scales::pretty_breaks(n_breaks)) +
		ggtitle("Ortho VITROS IgG") +
		xlab("Days") +
		ylab("log S/C value") -> plot_kinetics.VitrosIgG

## Plot sensitivity
as.data.frame(summary(fit_Stan.VitrosIgG)$summary) %>%
	tibble::rownames_to_column() %>%
	filter(str_detect(rowname, paste0("sensitivity_over_time\\["))) %>%
	mutate(severity_final_class = as.numeric(str_match_all(rowname, "(?<=\\[).+?(?=\\,)"))) %>%
	mutate(day = as.numeric(str_match_all(rowname, "(?<=\\,).+?(?=\\])"))) %>%
	mutate(severity_final_class = factor(severity_final_class, levels=c(1,2), labels=c("Not hospitalized", "Hospitalized"))) %>%
	# filter(day <= 140) %>%
	ggplot() +
		theme_bw() +
		theme(legend.position="none") +
		geom_vline(data=tmp_vertical, aes(xintercept=max_days), colour="purple") +
		geom_ribbon(aes(x=day, ymin=`2.5%`, ymax=`97.5%`, group=severity_final_class, fill=severity_final_class), alpha=0.2) +
		geom_line(aes(x=day, y=`50%`, group=severity_final_class, colour=severity_final_class), size=1) +
		scale_x_continuous(breaks=scales::pretty_breaks(n_breaks)) +
		scale_y_continuous(breaks=scales::pretty_breaks(5), limits=c(0,1)) +
		scale_fill_manual(values=c("Not hospitalized"="cornflowerblue", "Hospitalized"="red")) +
		scale_colour_manual(values=c("Not hospitalized"="cornflowerblue", "Hospitalized"="red")) +
		xlab("Days") +
		ylab("Sensitivity") +
		facet_wrap(~severity_final_class, scales="fixed", ncol=2) -> plot_Sens.VitrosIgG

## Putting it together
(plot_kinetics.Abbott + plot_kinetics.Roche + plot_kinetics.VitrosIgG) / (plot_Sens.Abbott + plot_Sens.Roche + plot_Sens.VitrosIgG)
