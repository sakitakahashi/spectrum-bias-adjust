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
library(fiftystater)

###################
## Load the data ##
###################

load("Results/Results_Italy_2groups.RData")
load("Results/Results_Italy_3groups.RData")

load("Results/Results_Spain_2groups.RData")
load("Results/Results_Spain_3groups.RData")

load("Results/Results_US_2groups.RData")

load("Results/Results_Manaus_2groups.RData")
load("Results/Results_Manaus_3groups.RData")

load("Results/Results_Japan_2groups_1assay.RData")
load("Results/Results_Japan_2groups_2assays.RData")

## For US: need to get Round 1 survey data as well
MAX_VALUE_US_2groups_4 <- MAX_VALUE_US_2groups
PR_spdf_2groups_4 <- PR_spdf_2groups
rm(MAX_VALUE_US_2groups, PR_spdf_2groups)

## Get round 1 data as well
dat_sero_US_track_2groups %>%
	filter(Round==1) -> survey_1_weighted_by_age_and_sex_2groups

## Re-link them back to states
dat_US_state_counts_2groups %>%
	filter(Round==1) %>%
	select(state, census_division) %>%
	left_join(survey_1_weighted_by_age_and_sex_2groups, by=c("census_division"="Census_Division")) -> survey_1_weighted_by_age_and_sex_final_2groups

survey_1_weighted_by_age_and_sex_final_2groups$state <- usdata::abbr2state(survey_1_weighted_by_age_and_sex_final_2groups$state)
survey_1_weighted_by_age_and_sex_final_2groups$state[49] <- "Puerto Rico"
survey_1_weighted_by_age_and_sex_final_2groups$state <- tolower(survey_1_weighted_by_age_and_sex_final_2groups$state)

## Add Puerto Rico
USA <- sf::st_read("Data/tl_2020_us_county.shp")
PR_counties <- USA %>% filter(STATEFP==72)
PR <- as(st_geometry(PR_counties), "Spatial")
PR_outline <- rgeos::gUnaryUnion(PR)
PR_spdf_2groups_1 <- broom::tidy(PR_outline)
PR_spdf_2groups_1$adj_sero_median <- survey_1_weighted_by_age_and_sex_final_2groups %>% filter(census_division=="South Atlantic & Puerto Rico") %>% slice(1) %>% select(adj_sero_median) %>% pull()

## What's the largest value to plot
survey_1_weighted_by_age_and_sex_final_2groups %>%
	select(-n_positive_raw, -n_tested_raw, -Round) %>%
	select_if(is.numeric) %>%
	max() + 0.01 -> MAX_VALUE_US_2groups_1

## Set plot limits
MAX_ITALY <- max(c(MAX_VALUE_ITALY_2groups, MAX_VALUE_ITALY_3groups))
MAX_US <- max(c(MAX_VALUE_US_2groups_1, MAX_VALUE_US_2groups_4))
MAX_OVERALL <- max(c(MAX_ITALY, MAX_VALUE_SPAIN_2groups, MAX_US))

###########
## Italy ##
###########

## Estimate single value for country
dat_sero_Italy_2groups %>%
	summarise(
		median = weighted.mean(adj_sero_median, Population_pct),
		lb = weighted.mean(adj_sero_lb, Population_pct),
		ub = weighted.mean(adj_sero_ub, Population_pct))

dat_sero_Italy_3groups %>%
	summarise(
		median = weighted.mean(adj_sero_median, Population_pct),
		lb = weighted.mean(adj_sero_lb, Population_pct),
		ub = weighted.mean(adj_sero_ub, Population_pct))

dat_sero_Italy_2groups %>%
	as.data.frame() %>%
	mutate(ratio = adj_sero_median/Seroprevalence, ratio_lb = adj_sero_lb/Seroprevalence, ratio_ub = adj_sero_ub/Seroprevalence) %>%
	arrange(desc(ratio))

Italy_df_2groups %>%
	mutate(Groups = "2 groups") %>%
	bind_rows(Italy_df_3groups %>% mutate(Groups = "3 groups (AS as separate)")) %>%
	ggplot() +
	geom_polygon(aes(x=long, y=lat, group=group, fill=adj_sero_median), color='black') +
	ggtitle("Italy") +
	coord_map() +
	theme_bw() +
	theme( ## Box only
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	 axis.ticks = element_blank(),
	axis.text=element_blank(),
	axis.line=element_blank()) +
	xlab("") + ylab("") +
	facet_wrap(.~Groups) +
	# scale_fill_gradientn("Median", limits=c(0,MAX_OVERALL), breaks=pretty_breaks(7), colours=brewer.pal(9, "YlOrRd")) +
	scale_fill_gradientn(name="Median",
		limits=c(0,MAX_OVERALL),
		values=c(0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, MAX_OVERALL)/MAX_OVERALL,
		colours=rev(brewer.pal(8, "RdYlBu"))) -> p_Italy_adj_median

###########
## Spain ##
###########

## Estimate single value for country
dat_sero_Spain_2groups %>%
	left_join(dat_demog_Spain, by=c("Provincia"="Province")) %>%
	group_by(Round) %>%
	summarise(
		median = weighted.mean(median, Population_pct),
		lb = weighted.mean(lb, Population_pct),
		ub = weighted.mean(ub, Population_pct),
		raw = weighted.mean(reported, Population_pct))

dat_sero_Spain_3groups %>%
	left_join(dat_demog_Spain, by=c("Provincia"="Province")) %>%
	group_by(Round) %>%
	summarise(
		median = weighted.mean(median, Population_pct),
		lb = weighted.mean(lb, Population_pct),
		ub = weighted.mean(ub, Population_pct),
		raw = weighted.mean(reported, Population_pct))

dat_sero_Spain_2groups %>%
	as.data.frame() %>%
	mutate(ratio = median/reported, ratio_lb = lb/reported, ratio_ub = ub/reported) %>%
	arrange(desc(ratio))

Spain_df_2groups %>%
	filter(lat >= 35) %>%
	select(-Seroprevalence_2, -adj_sero_2_median, -adj_sero_2_lb, -adj_sero_2_ub) %>%
	mutate(Round = "Round 1") %>%
	rename(Seroprevalence = Seroprevalence_1) %>%
	rename(adj_sero_median = adj_sero_1_median) %>%
	rename(adj_sero_lb = adj_sero_1_lb) %>%
	rename(adj_sero_ub = adj_sero_1_ub) -> tmp_Spain_1

Spain_df_2groups %>%
	filter(lat >= 35) %>%
	select(-Seroprevalence_1, -adj_sero_1_median, -adj_sero_1_lb, -adj_sero_1_ub) %>%
	mutate(Round = "Round 2") %>%
	rename(Seroprevalence = Seroprevalence_2) %>%
	rename(adj_sero_median = adj_sero_2_median) %>%
	rename(adj_sero_lb = adj_sero_2_lb) %>%
	rename(adj_sero_ub = adj_sero_2_ub) -> tmp_Spain_2

bind_rows(tmp_Spain_1, tmp_Spain_2) %>%
ggplot() +
	geom_polygon(aes(x=long, y=lat, group=group, fill=adj_sero_median), color='black') +
	ggtitle("Spain") +
	coord_map() +
	theme_bw() +
	theme( ## Box only
	# legend.position = "none",
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	 axis.ticks = element_blank(),
	axis.text=element_blank(),
	axis.line=element_blank()) +
	xlab("") + ylab("") +
	facet_wrap(.~Round) +
	# scale_fill_gradientn("Median", limits=c(0,MAX_OVERALL), breaks=pretty_breaks(7), colours=brewer.pal(9, "YlOrRd")) +
	scale_fill_gradientn(name="Median",
		limits=c(0,MAX_OVERALL),
		values=c(0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, MAX_OVERALL)/MAX_OVERALL,
		colours=rev(brewer.pal(8, "RdYlBu"))) -> p_Spain_adj_median_2groups

########
## US ##
########

ggplot(
	bind_rows(
		as.data.frame(survey_1_weighted_by_age_and_sex_final_2groups) %>% mutate(Round="Round 1"),
		as.data.frame(survey_4_weighted_by_age_and_sex_final_2groups) %>% mutate(Round="Round 4")), aes(map_id=state)) +
  geom_map(aes(fill=adj_sero_median), map = fifty_states, col="black", size=0.5) +
  coord_map() +
  geom_polygon(data=CD_1_2groups, aes(x=long, y=lat, group=group, map_id=1), fill=NA, colour="black", size=1, alpha=0) +  
  geom_polygon(data=CD_2_2groups, aes(x=long, y=lat, group=group, map_id=1), fill=NA, colour="black", size=1, alpha=0) +  
  geom_polygon(data=CD_3_2groups, aes(x=long, y=lat, group=group, map_id=1), fill=NA, colour="black", size=1, alpha=0) +  
  geom_polygon(data=CD_4_2groups, aes(x=long, y=lat, group=group, map_id=1), fill=NA, colour="black", size=1, alpha=0) +  
  geom_polygon(data=CD_5_2groups, aes(x=long, y=lat, group=group, map_id=1), fill=NA, colour="black", size=1, alpha=0) +  
  geom_polygon(data=CD_6_2groups, aes(x=long, y=lat, group=group, map_id=1), fill=NA, colour="black", size=1, alpha=0) +  
  geom_polygon(data=CD_7_2groups, aes(x=long, y=lat, group=group, map_id=1), fill=NA, colour="black", size=1, alpha=0) +  
  geom_polygon(data=CD_8_2groups, aes(x=long, y=lat, group=group, map_id=1), fill=NA, colour="black", size=1, alpha=0) +  
  geom_polygon(data=CD_9_2groups, aes(x=long, y=lat, group=group, map_id=1), fill=NA, colour="black", size=1, alpha=0) +  
  geom_polygon(data=PR_spdf_2groups_1 %>% mutate(Round="Round 1") %>% bind_rows(PR_spdf_2groups_4 %>% mutate(Round="Round 4")), aes(x=long-9+8, y=lat+7, group=group, map_id=1, fill=adj_sero_median), colour="black", size=0.5) +  
  geom_rect(aes(xmin=-123.9834, xmax=-112.3072, ymin=23.01948, ymax=31.24104), colour="black", fill=NA) +
  geom_rect(aes(xmin=-102.9840, xmax=-112.2745, ymin=23.01704, ymax=28.60362), colour="black", fill=NA) +
  geom_rect(aes(xmin=-78+8, xmax=-73+8, ymin=23, ymax=27), colour="black", fill=NA) +
  expand_limits(x = fifty_states$long, y = fifty_states$lat) +
	ggtitle("United States") +
	theme_bw() +
	theme( ## Box only
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	 axis.ticks = element_blank(),
	axis.text=element_blank(),
	axis.line=element_blank()) +
	xlab("") + ylab("") +
	facet_wrap(.~Round) +
	# scale_fill_gradientn("Median", limits=c(0,MAX_OVERALL), breaks=pretty_breaks(7), colours=brewer.pal(9, "YlOrRd")) +
	scale_fill_gradientn(name="Median",
		limits=c(0,MAX_OVERALL),
		values=c(0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, MAX_OVERALL)/MAX_OVERALL,
		colours=rev(brewer.pal(8, "RdYlBu"))) -> p_US_adj_median_2groups

############
## Manaus ##
############

dat_sero_Manaus_AGG_2groups %>%
	mutate(Groups = "2 groups") %>%
	bind_rows(dat_sero_Manaus_AGG_3groups %>% mutate(Groups = "3 groups\n(AS as separate)")) %>%
	mutate(Groups = factor(Groups, levels=c("3 groups\n(AS as separate)", "2 groups"))) %>%
	ggplot() +
		theme_bw() +
		# theme(legend.position="none") +
		geom_linerange(aes(x=clean_date-days(14), y=adj_sero_median, ymin=adj_sero_lb, ymax=adj_sero_ub, group=Groups), size=0.5, position = position_dodge(width=5)) +
		geom_line(aes(x=clean_date-days(14), y=adj_sero_median, colour=Groups), position = position_dodge(width=5)) +
		geom_point(aes(x=clean_date-days(14), y=adj_sero_median, fill=Groups), position = position_dodge(width=5), size=3, shape=21) +
		## Add in raw data
		geom_point(aes(x=clean_date-days(14), y=raw_prev), position=position_dodge(width=5), size=2, shape=23, fill="black") +
		scale_x_date(date_labels="%b %Y", date_breaks="1 month") +
		xlab("") +
		ggtitle("Manaus, Brazil") +
		theme(
			axis.text.x=element_text(angle=45, hjust=1),
			panel.grid.minor.x = element_blank()) +
		ylim(0,1) +
		scale_fill_manual(values=c("2 groups"="deeppink3", "3 groups\n(AS as separate)"="aquamarine3")) +
		scale_colour_manual(values=c("2 groups"="deeppink3", "3 groups\n(AS as separate)"="aquamarine3")) +
		ylab("Cumulative incidence") -> p_Manaus

dat_sero_Manaus_AGG_2groups %>%
	mutate(ratio = adj_sero_median/raw_prev, ratio_lb = adj_sero_lb/raw_prev, ratio_ub = adj_sero_ub/raw_prev) %>%
	arrange(desc(ratio))

###########
## Japan ##
###########

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
	left_join(tmp_raw %>% rename(raw=value) %>% select(Prefecture, Assay, raw), by=c("Prefecture", "Assay")) %>%
	mutate(Assay = factor(Assay, levels=c("Abbott","Roche","Bivariate"), labels=c("Abbott\nonly","Roche\nonly","Both\nassays"))) %>%
	mutate(Prefecture = factor(Prefecture, levels=c("Tokyo","Osaka","Miyagi","Aichi","Fukuoka"))) %>%
	ggplot() +
		geom_linerange(aes(x=Assay, y=adj_median, ymin=adj_lb, ymax=adj_ub), size=0.5) +
		geom_point(aes(x=Assay, y=adj_median, fill="Adjusted"), shape=21, size=3, stroke=0.5) +
		## Add in raw data
		geom_point(aes(x=Assay, y=raw, fill="Raw"), size=2, shape=23) +
		facet_wrap(.~Prefecture, scales="fixed", ncol=3) +
		theme_bw() +
		xlab("") +
		scale_fill_manual("Value", values=c("Raw"="black", "Adjusted"="deeppink3")) +
		ylab("Cumulative incidence") +
		ggtitle("Japan") -> p_Japan

## Put it all together
(p_Italy_adj_median + p_Spain_adj_median_2groups) + plot_layout(guides = 'collect')
p_US_adj_median_2groups + plot_layout(guides = 'collect')
p_Manaus + p_Japan
