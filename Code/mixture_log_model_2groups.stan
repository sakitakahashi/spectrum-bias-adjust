data {
	
	int N;									// Number of observations
	int J;									// Number of individuals
	
	int patient_ID_index[N];
	
	real<lower=0.0> time_since[N];
	
	real log_SC_value[N];
	
	int severity_by_obs[N];
	
	int patient_ID_index_by_severity[N];
	
	int J_1;								// Number of individuals in severity group 1: non-hosp
	int J_2;								// Number of individuals in severity group 2 hosp
	
	int N_rep;
	
	real CUTOFF;
	
	int DAYS_FOR_SENS;
	
	int N_UNIQUE_LAMBDA;					// Either 1 or 2
	
	int severity_by_ind[J];
	int patient_ID_index_by_severity_by_ind[J];
	
}
parameters {
	
	real lambda[N_UNIQUE_LAMBDA];			// Order: low, high
	
	ordered[2] log_beta_severity;			// Order: low, high
	
	vector[J_1] log_intercept_raw_1;
	
	vector[J_2] log_intercept_raw_2;
	
	real<lower=0.0> sigma;
	
	positive_ordered[2] tau_backwards;		// Order: high, low
	
	real<lower=0.0, upper=1.0> theta[2]; 	// Probability of being a low (theta[s]) or high (1-theta[s]) responder in group s
	
}
transformed parameters {
	
	real<lower=0.0> tau[2];					// Order: low, high
	
	tau[1] = tau_backwards[2];
	tau[2] = tau_backwards[1];
	
}
model {
	
	// Likelihood part 1: intercepts for each group
	for (j in 1:J_1) {
			
			target += log_mix(theta[1],
				normal_lpdf(log_intercept_raw_1[j] | log_beta_severity[1], tau[1]),
				normal_lpdf(log_intercept_raw_1[j] | log_beta_severity[2], tau[2]));
			
	}
	
	for (j in 1:J_2) {
			
			target += log_mix(theta[2],
				normal_lpdf(log_intercept_raw_2[j] | log_beta_severity[1], tau[1]),
				normal_lpdf(log_intercept_raw_2[j] | log_beta_severity[2], tau[2]));
			
	}
	
	// Likelihood part 2: data
	for (n in 1:N) {
		
		if(severity_by_obs[n]==1) {
			
			if(N_UNIQUE_LAMBDA==1) target += normal_lpdf(log_SC_value[n] | log_intercept_raw_1[patient_ID_index_by_severity[n]] - (lambda[1] * time_since[n]), sigma);
			
			else if(N_UNIQUE_LAMBDA==2) {
				
				target += log_mix(theta[1],
					normal_lpdf(log_SC_value[n] | log_intercept_raw_1[patient_ID_index_by_severity[n]] - (lambda[1] * time_since[n]), sigma),
					normal_lpdf(log_SC_value[n] | log_intercept_raw_1[patient_ID_index_by_severity[n]] - (lambda[2] * time_since[n]), sigma));
				
			}
			
		}
		
		if(severity_by_obs[n]==2) {
			
			if(N_UNIQUE_LAMBDA==1) target += normal_lpdf(log_SC_value[n] | log_intercept_raw_2[patient_ID_index_by_severity[n]] - (lambda[1] * time_since[n]), sigma);
			
			else if(N_UNIQUE_LAMBDA==2) {
				
				target += log_mix(theta[2],
					normal_lpdf(log_SC_value[n] | log_intercept_raw_2[patient_ID_index_by_severity[n]] - (lambda[1] * time_since[n]), sigma),
					normal_lpdf(log_SC_value[n] | log_intercept_raw_2[patient_ID_index_by_severity[n]] - (lambda[2] * time_since[n]), sigma));
				
			}
			
		}
		
	}
	
	// Priors
	if(N_UNIQUE_LAMBDA==1) target += normal_lpdf(lambda[1] | 0.0, 0.01);
	if(N_UNIQUE_LAMBDA==2) {
		
		target += normal_lpdf(lambda[1] | 0.0, 0.01);
		target += normal_lpdf(lambda[2] | 0.0, 0.01);
		
	}
	
	target += normal_lpdf(log_beta_severity[1] | 0.0, 2.5);
	target += normal_lpdf(log_beta_severity[2] | 0.0, 2.5);
	
	target += normal_lpdf(tau_backwards[1] | 0.0, 1.0);
	target += normal_lpdf(tau_backwards[2] | 0.0, 1.0);
	
	target += normal_lpdf(theta[1] | 0.5, 0.5);
	target += normal_lpdf(theta[2] | 0.5, 0.5);
	
}
generated quantities {
	
	matrix<lower=0.0, upper=1.0>[2,(DAYS_FOR_SENS+1)] sensitivity_over_time;	// Sensitivity: 2 classes of severity (non-hosp, hosp)
	
	real log_lik[N];
	
	matrix[J,(DAYS_FOR_SENS+1)] log_SC_value_fitted;
	
	// Calculate time-varying sensitivity by symptoms
	for(k in 0:DAYS_FOR_SENS) {
		
		vector[N_rep] which_group_1_low;
		vector[N_rep] which_group_2_low;
		
		vector[N_rep] intercept_new_1_tau;
		vector[N_rep] intercept_new_2_tau;
		
		vector[N_rep] intercept_new_sigma;
		
		vector[N_rep] Y_tilde_1;
		vector[N_rep] Y_tilde_2;
		
		vector[N_rep] Y_tilde_binary_1;
		vector[N_rep] Y_tilde_binary_2;
		
		// Simulate new data
		for(x in 1:N_rep) {
			
			// Assign the observation for each class as a low or high responder
			which_group_1_low[x] = bernoulli_rng(theta[1]);
			which_group_2_low[x] = bernoulli_rng(theta[2]);
			
			// Generate other random variables
			intercept_new_1_tau[x] = normal_rng(0.0, tau[1]);
			intercept_new_2_tau[x] = normal_rng(0.0, tau[2]);
			
			intercept_new_sigma[x] = normal_rng(0.0, sigma);
			
			// Generate continuous Y
			if(N_UNIQUE_LAMBDA==1) {
				
				Y_tilde_1[x] = which_group_1_low[x] == 1 ?
					log_beta_severity[1] + intercept_new_1_tau[x] + intercept_new_sigma[x] - (lambda[1] * (k * 1.0)) :
					log_beta_severity[2] + intercept_new_2_tau[x] + intercept_new_sigma[x] - (lambda[1] * (k * 1.0));
				
				Y_tilde_2[x] = which_group_2_low[x] == 1 ?
					log_beta_severity[1] + intercept_new_1_tau[x] + intercept_new_sigma[x] - (lambda[1] * (k * 1.0)) :
					log_beta_severity[2] + intercept_new_2_tau[x] + intercept_new_sigma[x] - (lambda[1] * (k * 1.0));
				
			}
			
			else if(N_UNIQUE_LAMBDA==2) {
				
				Y_tilde_1[x] = which_group_1_low[x] == 1 ?
					log_beta_severity[1] + intercept_new_1_tau[x] + intercept_new_sigma[x] - (lambda[1] * (k * 1.0)) :
					log_beta_severity[2] + intercept_new_2_tau[x] + intercept_new_sigma[x] - (lambda[2] * (k * 1.0));
				
				Y_tilde_2[x] = which_group_2_low[x] == 1 ?
					log_beta_severity[1] + intercept_new_1_tau[x] + intercept_new_sigma[x] - (lambda[1] * (k * 1.0)) :
					log_beta_severity[2] + intercept_new_2_tau[x] + intercept_new_sigma[x] - (lambda[2] * (k * 1.0));
				
			}
			
			// Transform Y to binary results based on cutoff
			Y_tilde_binary_1[x] = Y_tilde_1[x] >= CUTOFF ? 1.0 : 0.0;
			Y_tilde_binary_2[x] = Y_tilde_2[x] >= CUTOFF ? 1.0 : 0.0;
			
		}
		
		// Calculate sensitivity
		sensitivity_over_time[1,k+1] = sum(Y_tilde_binary_1) / (1.0 * N_rep);
		sensitivity_over_time[2,k+1] = sum(Y_tilde_binary_2) / (1.0 * N_rep);
		
	}
	
	// Calculate pointwise likelihood
	for (n in 1:N) {
		
		if(severity_by_obs[n]==1) {
			
			if(N_UNIQUE_LAMBDA==1) log_lik[n] = normal_lpdf(log_SC_value[n] | log_intercept_raw_1[patient_ID_index_by_severity[n]] - (lambda[1] * time_since[n]), sigma);
			
			else if(N_UNIQUE_LAMBDA==2) {
				
				log_lik[n] = log_mix(theta[1],
					normal_lpdf(log_SC_value[n] | log_intercept_raw_1[patient_ID_index_by_severity[n]] - (lambda[1] * time_since[n]), sigma),
					normal_lpdf(log_SC_value[n] | log_intercept_raw_1[patient_ID_index_by_severity[n]] - (lambda[2] * time_since[n]), sigma));
				
			}
			
		}
		
		if(severity_by_obs[n]==2) {
			
			if(N_UNIQUE_LAMBDA==1) log_lik[n] = normal_lpdf(log_SC_value[n] | log_intercept_raw_2[patient_ID_index_by_severity[n]] - (lambda[1] * time_since[n]), sigma);
			
			else if(N_UNIQUE_LAMBDA==2) {
				
				log_lik[n] = log_mix(theta[2],
					normal_lpdf(log_SC_value[n] | log_intercept_raw_2[patient_ID_index_by_severity[n]] - (lambda[1] * time_since[n]), sigma),
					normal_lpdf(log_SC_value[n] | log_intercept_raw_2[patient_ID_index_by_severity[n]] - (lambda[2] * time_since[n]), sigma));
				
			}
			
		}
		
	}
	
	// The posterior predictive distribution for each person by time
	for(j in 1:J) {
		
		for(k in 0:DAYS_FOR_SENS) {
			
			if(severity_by_ind[j]==1) log_SC_value_fitted[j,k+1] = normal_rng(log_intercept_raw_1[patient_ID_index_by_severity_by_ind[j]] - (lambda[1] * (k * 1.0)), sigma);
			if(severity_by_ind[j]==2) log_SC_value_fitted[j,k+1] = normal_rng(log_intercept_raw_2[patient_ID_index_by_severity_by_ind[j]] - (lambda[1] * (k * 1.0)), sigma);
			
		}
		
	}
	
}
