## This leverages standalone generated quantities:
## https://mc-stan.org/rstan/reference/stanmodel-method-gqs.html

#######################################################################
## Function to calculate weighted Se by time, adjusting for severity ##
#######################################################################

## Version 8 of the kinetics code estimates time-varying Se up to 365 days
## Version 7 of the kinetics code estimates time-varying Se up to 200 days

calc_wtd_Se_by_time_adjusting_for_severity <-'
data {
	
	real prop[3];
	
}
parameters {
	
	matrix<lower=0.0, upper=1.0>[3,366] sensitivity_over_time;
	
}
generated quantities {
	
	real sensitivity_wtd[366];
	
	for(i in 1:366) sensitivity_wtd[i] = (prop[1]*sensitivity_over_time[1,i]) + (prop[2]*sensitivity_over_time[2,i]) + (prop[3]*sensitivity_over_time[3,i]);
	
}
'

calc_wtd_Se_by_time_adjusting_for_severity_2groups <-'
data {
	
	real prop[2];
	
}
parameters {
	
	matrix<lower=0.0, upper=1.0>[2,366] sensitivity_over_time;
	
}
generated quantities {
	
	real sensitivity_wtd[366];
	
	for(i in 1:366) sensitivity_wtd[i] = (prop[1]*sensitivity_over_time[1,i]) + (prop[2]*sensitivity_over_time[2,i]);
	
}
'


#####################################################################################
## Function to calculate weighted Se by time, adjusting for severity & time series ##
#####################################################################################

calc_wtd_Se_by_time_adjusting_for_severity_and_ts <-'
data {
	
	int N;
	real prop_infections[N];
	
}
parameters {
	
	real sensitivity_wtd[N];
	
}
generated quantities {
	
	real sensitivity_wtd_with_ts;
	
	sensitivity_wtd_with_ts = 0.0;
	
	for(i in 1:N) sensitivity_wtd_with_ts += sensitivity_wtd[i] * prop_infections[N+1-i];
	
}
'

calc_wtd_Se_by_time_adjusting_for_severity_and_ts_3_assays <-'
data {
	
	int N;
	real prop_infections[N];
	
}
parameters {
	
	real sensitivity_wtd_Abbott[N];
	real sensitivity_wtd_VitrosIgG[N];
	real sensitivity_wtd_Roche[N];
	
}
generated quantities {
	
	real sensitivity_wtd_Abbott_with_ts;
	real sensitivity_wtd_VitrosIgG_with_ts;
	real sensitivity_wtd_Roche_with_ts;
	
	sensitivity_wtd_Abbott_with_ts = 0.0;
	sensitivity_wtd_VitrosIgG_with_ts = 0.0;
	sensitivity_wtd_Roche_with_ts = 0.0;
	
	for(i in 1:N) {
		
		sensitivity_wtd_Abbott_with_ts += sensitivity_wtd_Abbott[i] * prop_infections[N+1-i];
		sensitivity_wtd_VitrosIgG_with_ts += sensitivity_wtd_VitrosIgG[i] * prop_infections[N+1-i];
		sensitivity_wtd_Roche_with_ts += sensitivity_wtd_Roche[i] * prop_infections[N+1-i];
		
	}
}
'

#############################################################################################################################
## Function to calculate weighted Se by time (adjusting for severity & time series) and adjusted sero-prevalence using R-G ##
#############################################################################################################################

calc_wtd_Se_by_time_adjusting_for_severity_and_ts_with_seroprev_RG <-'
data {
int N;
real prop_infections[N];
real prev_raw;
real Sp;
}
parameters {
real sensitivity_wtd[N];
}
generated quantities {
real sensitivity_wtd_with_ts;
real prev_adj;

sensitivity_wtd_with_ts = 0.0;
for(i in 1:N) sensitivity_wtd_with_ts += sensitivity_wtd[i] * prop_infections[N+1-i];

prev_adj = (prev_raw + Sp - 1.0) / (sensitivity_wtd_with_ts + Sp - 1.0);

}
'

#############################################################
## Function to obtain binomial estimates of seroprevalence ##
#############################################################

calc_seroprev_binomial <- '
data {
	
	int y_sample;
	int n_sample;
	
	real mu_Se;
	real sigma_Se;
	
	real Sp;
	
}
parameters {
	
	real<lower=0.0, upper=1.0> prev_adj;
	real<lower=0.0, upper=1.0> Se;
	
}
model {
	
	real p_sample = prev_adj*Se + (1.0-prev_adj)*(1-Sp);
	
	y_sample ~ binomial(n_sample, p_sample);
	
	Se ~ normal(mu_Se, sigma_Se);
	
}'

## Alternatively
calc_seroprev_binomial_logit_scale <- '
data {
	
	int y_sample;
	int n_sample;
	
	real mu_Se;
	real sigma_Se;
	
	real Sp;
	
}
parameters {
	
	real logit_prev_adj;
	real<lower=0.0, upper=1.0> Se;
	
}
transformed parameters {
	
	real<lower=0.0, upper=1.0> prev_adj = inv_logit(logit_prev_adj);
	
}
model {
	
	real p_sample = prev_adj*Se + (1.0-prev_adj)*(1-Sp);
	
	y_sample ~ binomial(n_sample, p_sample);
	
	Se ~ normal(mu_Se, sigma_Se);
	
}'

calc_seroprev_binomial_3_assays <-'
data {
	
	int y_sample_Abbott;
	int n_sample_Abbott;
	
	int y_sample_VitrosIgG;
	int n_sample_VitrosIgG;
	
	int y_sample_Roche;
	int n_sample_Roche;
	
	real mu_Se_Abbott;
	real sigma_Se_Abbott;
	
	real mu_Se_VitrosIgG;
	real sigma_Se_VitrosIgG;
	
	real mu_Se_Roche;
	real sigma_Se_Roche;
	
	//real prop_test[3];
	
	real Sp_Abbott;
	real Sp_VitrosIgG;
	real Sp_Roche;
	
}
parameters {
	
	real<lower=0.0, upper=1.0> prev_adj;
	
	real<lower=0.0, upper=1.0> Se_Abbott;
	real<lower=0.0, upper=1.0> Se_VitrosIgG;
	real<lower=0.0, upper=1.0> Se_Roche;
	
}
transformed parameters {
	
	real<lower=0.0, upper=1.0> p_sample[3];
	
	p_sample[1] = prev_adj*Se_Abbott + (1.0-prev_adj)*(1.0-Sp_Abbott);
	p_sample[2] = prev_adj*Se_VitrosIgG + (1.0-prev_adj)*(1.0-Sp_VitrosIgG);
	p_sample[3] = prev_adj*Se_Roche + (1.0-prev_adj)*(1.0-Sp_Roche);
	
}
model{
	
	y_sample_Abbott ~ binomial(n_sample_Abbott, p_sample[1]);
	y_sample_VitrosIgG ~ binomial(n_sample_VitrosIgG, p_sample[2]);
	y_sample_Roche ~ binomial(n_sample_Roche, p_sample[3]);
	
	Se_Abbott ~ normal(mu_Se_Abbott, sigma_Se_Abbott);
	Se_VitrosIgG ~ normal(mu_Se_VitrosIgG, sigma_Se_VitrosIgG);
	Se_Roche ~ normal(mu_Se_Roche, sigma_Se_Roche);
	
}
'

calc_seroprev_binomial_adjusting_over_n_groups <- '
data {
	
	int n_groups;
	vector[n_groups] weight;
	
	int y_sample[n_groups];
	int n_sample[n_groups];
	
	real mu_Se[n_groups];
	real sigma_Se[n_groups];
	
	real Sp;
	
}
parameters {
	
	vector<lower=0.0, upper=1.0>[n_groups] prev_adj;
	vector<lower=0.0, upper=1.0>[n_groups] Se;
	
}
model {
	
	real p_sample[n_groups];
	
	for(g in 1:n_groups) {
	
		p_sample[g] = prev_adj[g]*Se[g] + (1.0-prev_adj[g])*(1-Sp);
		
		y_sample[g] ~ binomial(n_sample[g], p_sample[g]);
		
		Se[g] ~ normal(mu_Se[g], sigma_Se[g]);
		
	}
	
}
generated quantities {
	
	real prev_adj_overall;
	
	prev_adj_overall = dot_product(weight, prev_adj);
	
}
'

calc_seroprev_binomial_3_assays_adjusting_over_n_groups_per_assay <-'
data {
	
	int n_groups;				// 4 age groups, 2 sex = 8 groups total
	vector[n_groups] weight;	// 4 age groups, 2 sex = 8 groups total
	
	int y_sample_Abbott[n_groups];
	int n_sample_Abbott[n_groups];
	
	int y_sample_VitrosIgG[n_groups];
	int n_sample_VitrosIgG[n_groups];
	
	int y_sample_Roche[n_groups];
	int n_sample_Roche[n_groups];
	
	real mu_Se_Abbott[n_groups];
	real sigma_Se_Abbott[n_groups];
	
	real mu_Se_VitrosIgG[n_groups];
	real sigma_Se_VitrosIgG[n_groups];
	
	real mu_Se_Roche[n_groups];
	real sigma_Se_Roche[n_groups];
	
	//real prop_test[3];
	
	real Sp_Abbott;
	real Sp_VitrosIgG;
	real Sp_Roche;
	
}
parameters {
	
	vector<lower=0.0, upper=1.0>[n_groups] prev_adj;
	
	vector<lower=0.0, upper=1.0>[n_groups] Se_Abbott;
	vector<lower=0.0, upper=1.0>[n_groups] Se_VitrosIgG;
	vector<lower=0.0, upper=1.0>[n_groups] Se_Roche;
	
}
transformed parameters {
	
	real<lower=0.0, upper=1.0> p_sample[3,n_groups];
	
	for(g in 1:n_groups) {
		
		p_sample[1,g] = prev_adj[g]*Se_Abbott[g] + (1.0-prev_adj[g])*(1.0-Sp_Abbott);
		p_sample[2,g] = prev_adj[g]*Se_VitrosIgG[g] + (1.0-prev_adj[g])*(1.0-Sp_VitrosIgG);
		p_sample[3,g] = prev_adj[g]*Se_Roche[g] + (1.0-prev_adj[g])*(1.0-Sp_Roche);
		
	}
	
}
model{
	
	for(g in 1:n_groups) {
		
		y_sample_Abbott[g] ~ binomial(n_sample_Abbott[g], p_sample[1,g]);
		y_sample_VitrosIgG[g] ~ binomial(n_sample_VitrosIgG[g], p_sample[2,g]);
		y_sample_Roche[g] ~ binomial(n_sample_Roche[g], p_sample[3,g]);
		
		Se_Abbott[g] ~ normal(mu_Se_Abbott[g], sigma_Se_Abbott[g]);
		Se_VitrosIgG[g] ~ normal(mu_Se_VitrosIgG[g], sigma_Se_VitrosIgG[g]);
		Se_Roche[g] ~ normal(mu_Se_Roche[g], sigma_Se_Roche[g]);
		
	}
	
}
generated quantities {
	
	real prev_adj_overall;
	
	prev_adj_overall = dot_product(weight, prev_adj);
	
}
'
