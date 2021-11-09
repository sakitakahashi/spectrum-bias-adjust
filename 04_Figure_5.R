library(readstata13)
library(dplyr)
library(data.table)
library(ggplot2)
library(tibble)
library(knitr)
library(tidyr)
library(cowplot)
library(stringr)
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

## Italy
load("Results/Results_Italy_2groups.RData")

dat_sero_Italy_2groups %>%
	ggplot(aes(x=Seroprevalence, y=adj_sero_median/Seroprevalence)) +
		geom_point(position=position_dodge2(width=0.0005), size=1.5) +
		geom_linerange(aes(ymin=adj_sero_lb/Seroprevalence, ymax=adj_sero_ub/Seroprevalence), position=position_dodge2(width=0.0005)) +
		geom_hline(yintercept=1, col="black", linetype="dashed") +
		theme_bw() +
		xlab("Raw seroprevalence") +
		ylab("Ratio, cumulative incidence to raw seroprevalence") +
		ggtitle("Italy") -> plot.Italy

## Spain
load("Results/Results_Spain_2groups.RData")

dat_sero_Spain_2groups %>%
	ggplot(aes(x=reported, y=median/reported)) +
		geom_point(position=position_dodge2(width=0.0005), size=1.5) +
		geom_linerange(aes(ymin=lb/reported, ymax=ub/reported), position=position_dodge2(width=0.0005)) +
		geom_hline(yintercept=1, col="black", linetype="dashed") +
		theme_bw() +
		xlab("Raw seroprevalence") +
		ylab("Ratio, cumulative incidence to raw seroprevalence") +
		ggtitle("Spain") +
		facet_wrap(.~Round, ncol=1) -> plot.Spain

## US
load("Results/Results_US_2groups.RData")

dat_sero_US_track_2groups %>%
	mutate(Round = factor(Round)) %>%
	rename(`Census Division` = Census_Division) %>%
	ggplot(aes(x=seroprevalence_raw, y=adj_sero_median/seroprevalence_raw)) +
		geom_linerange(aes(ymin=adj_sero_lb/seroprevalence_raw, ymax=adj_sero_ub/seroprevalence_raw), position=position_dodge2(width=0.0005), size=0.5) +
		geom_point(aes(fill=Round), position=position_dodge2(width=0.0005), shape=21, size=3) +
		geom_hline(yintercept=1, col="black", linetype="dashed") +
		theme_bw() +
		xlab("Raw seroprevalence") +
		ylab("Ratio, cumulative incidence to raw seroprevalence") +
		ggtitle("United States") +
		facet_wrap(.~`Census Division`, ncol=3, scales="free_x") +
		scale_fill_brewer(palette="Spectral") +
		theme(panel.spacing.x = unit(7, "mm")) -> plot.US

## Manaus
load("Results/Results_Manaus_2groups.RData")

dat_sero_Manaus_AGG_2groups %>%
	ggplot(aes(x=raw_prev, y=adj_sero_median/raw_prev)) +
		geom_linerange(aes(ymin=adj_sero_lb/raw_prev, ymax=adj_sero_ub/raw_prev), position=position_dodge2(width=0.0005), size=0.5) +
		geom_point(aes(fill=as.factor(clean_date - days(14))), position=position_dodge2(width=0.0005), shape=23, size=3) +
		geom_hline(yintercept=1, col="black", linetype="dashed") +
		theme_bw() +
		xlab("Raw seroprevalence") +
		ylab("Ratio, cumulative incidence to raw seroprevalence") +
		ggtitle("Manaus, Brazil") +
		theme(legend.position = c(0.8, 0.8)) +
		guides(fill=guide_legend(title="Month", ncol=2)) +
		scale_fill_brewer("Month", palette="Spectral", labels=as.Date_origin) -> plot.Manaus

## Japan
load("Results/Results_Japan_2groups_1assay.RData")
load("Results/Results_Japan_2groups_2assays.RData")

dat_sero_Japan_univariate_2groups %>%
	pivot_longer(cols=Abbott_raw:Roche_ub) %>%
	separate(name, sep="_", into=c("Assay","metric")) %>%
	filter(metric=="raw") -> tmp_raw_uni

dat_sero_Japan_bivariate_2groups %>%
	pivot_longer(cols=bivariate_raw:bivariate_ub) %>%
	separate(name, sep="_", into=c("Assay","metric")) %>%
	filter(metric=="raw") %>%
	mutate(metric = "raw") %>%
	mutate(Assay = "Bivariate") -> tmp_raw_bi

bind_rows(tmp_raw_uni, tmp_raw_bi) -> tmp_raw

dat_sero_Japan_univariate_2groups %>%
	mutate(Assay="Abbott") %>%
	mutate(metric="adjusted") %>%
	select(Prefecture, Assay, Abbott_median, Abbott_lb, Abbott_ub, metric) %>%
	rename(adj_median = Abbott_median) %>%
	rename(adj_lb = Abbott_lb) %>%
	rename(adj_ub = Abbott_ub) -> tmp_adj_uni_Abbott

dat_sero_Japan_univariate_2groups %>%
	mutate(Assay="Roche") %>%
	mutate(metric="adjusted") %>%
	select(Prefecture, Assay, Roche_median, Roche_lb, Roche_ub, metric) %>%
	rename(adj_median = Roche_median) %>%
	rename(adj_lb = Roche_lb) %>%
	rename(adj_ub = Roche_ub) -> tmp_adj_uni_Roche

dat_sero_Japan_bivariate_2groups %>%
	mutate(Assay="Bivariate") %>%
	mutate(metric="adjusted") %>%
	select(Prefecture, Assay, bivariate_median, bivariate_lb, bivariate_ub, metric) %>%
	rename(adj_median = bivariate_median) %>%
	rename(adj_lb = bivariate_lb) %>%
	rename(adj_ub = bivariate_ub) -> tmp_adj_bi

bind_rows(tmp_adj_uni_Abbott, tmp_adj_uni_Roche, tmp_adj_bi) -> tmp_adj

tmp_adj %>%
	select(-metric) %>%
	left_join(tmp_raw %>% rename(raw=value) %>% select(Prefecture, Assay, raw), by=c("Prefecture", "Assay")) -> dat_to_plot_Japan_2groups

dat_to_plot_Japan_2groups %>%
	filter(Assay == "Bivariate") %>%
	droplevels() %>%
	ggplot(aes(x=raw, y=adj_median/raw)) +
		geom_linerange(aes(ymin=adj_lb/raw, ymax=adj_ub/raw), position=position_dodge2(width=0.0005), size=0.5) +
		geom_point(aes(fill=Prefecture), position=position_dodge2(width=0.0005), shape=22, size=3) +
		geom_hline(yintercept=1, col="black", linetype="dashed") +
		theme_bw() +
		xlab("Raw seroprevalence") +
		ylab("Ratio, cumulative incidence to raw seroprevalence") +
		ggtitle("Japan") +
		theme(legend.position = c(0.8, 0.8)) +
		scale_fill_brewer(palette="Spectral") -> plot.Japan

## Put it all together
(plot.Italy + plot.Spain + plot.Manaus + plot.Japan | plot.US) + plot_annotation(tag_levels = 'A')
