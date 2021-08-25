data {
	
	int N;											// Number of observations
	int J;											// Number of individuals
	
	int patient_ID_index[N];
	
	real<lower=0.0> time_since[N];
	
	matrix[N,2] log_SC_value;
	
	int severity_by_obs[N];
	
	int patient_ID_index_by_severity[N];
	
	int J_1;								// Number of individuals in severity group 1: AS
	int J_2;								// Number of individuals in severity group 2: S & non-hosp
	int J_3;								// Number of individuals in severity group 3: hosp
	
	int N_rep;
	
	real CUTOFF[2];
	
	int DAYS_FOR_SENS;
	
}
parameters {
	
	real lambda_Abbott;								// Order: low, high [assuming slopes are shared between low and high]
	real lambda_Roche;								// Order: low, high [assuming slopes are shared between low and high]
	
	ordered[2] log_beta_severity_Abbott;			// Order: low, high
	ordered[2] log_beta_severity_Roche;				// Order: low, high
	
	vector[J_1] log_intercept_raw_1_Abbott;
	vector[J_1] log_intercept_raw_1_Roche;
	
	vector[J_2] log_intercept_raw_2_Abbott;
	vector[J_2] log_intercept_raw_2_Roche;
	
	vector[J_3] log_intercept_raw_3_Abbott;
	vector[J_3] log_intercept_raw_3_Roche;
	
	cov_matrix[2] Sigma;

	positive_ordered[2] tau_backwards_Abbott;		// Order: high, low
	positive_ordered[2] tau_backwards_Roche;		// Order: high, low
	
	real<lower=0.0, upper=1.0> theta_1_Abbott; 		// Probability of being a low (theta) or high (1-theta) responder in group 1
	real<lower=0.0, upper=1.0> theta_1_Roche; 		// Probability of being a low (theta) or high (1-theta) responder in group 1
	
	real<lower=0.0, upper=1.0> theta_2_Abbott; 		// Probability of being a low (theta) or high (1-theta) responder in group 2
	real<lower=0.0, upper=1.0> theta_2_Roche; 		// Probability of being a low (theta) or high (1-theta) responder in group 2
	
	real<lower=0.0, upper=1.0> theta_3_Abbott; 		// Probability of being a low (theta) or high (1-theta) responder in group 3
	real<lower=0.0, upper=1.0> theta_3_Roche; 		// Probability of being a low (theta) or high (1-theta) responder in group 3
	
}
transformed parameters {
	
	real<lower=0.0> tau_Abbott[2];					// Order: low, high
	real<lower=0.0> tau_Roche[2];					// Order: low, high
	
	tau_Abbott[1] = tau_backwards_Abbott[2];
	tau_Abbott[2] = tau_backwards_Abbott[1];
	
	tau_Roche[1] = tau_backwards_Roche[2];
	tau_Roche[2] = tau_backwards_Roche[1];
	
}
model {
	
	// Likelihood part 1: intercepts for groups 1 & 2
	for (j in 1:J_1) {
		
		target += log_mix(theta_1_Abbott,
			normal_lpdf(log_intercept_raw_1_Abbott[j] | log_beta_severity_Abbott[1], tau_Abbott[1]),
			normal_lpdf(log_intercept_raw_1_Abbott[j] | log_beta_severity_Abbott[2], tau_Abbott[2]));
		
		target += log_mix(theta_1_Roche,
			normal_lpdf(log_intercept_raw_1_Roche[j] | log_beta_severity_Roche[1], tau_Roche[1]),
			normal_lpdf(log_intercept_raw_1_Roche[j] | log_beta_severity_Roche[2], tau_Roche[2]));
		
	}
	
	for (j in 1:J_2) {
		
		target += log_mix(theta_2_Abbott,
			normal_lpdf(log_intercept_raw_2_Abbott[j] | log_beta_severity_Abbott[1], tau_Abbott[1]),
			normal_lpdf(log_intercept_raw_2_Abbott[j] | log_beta_severity_Abbott[2], tau_Abbott[2]));
		
		target += log_mix(theta_2_Roche,
			normal_lpdf(log_intercept_raw_2_Roche[j] | log_beta_severity_Roche[1], tau_Roche[1]),
			normal_lpdf(log_intercept_raw_2_Roche[j] | log_beta_severity_Roche[2], tau_Roche[2]));
		
	}
	
	for (j in 1:J_3) {
		
		target += log_mix(theta_3_Abbott,
			normal_lpdf(log_intercept_raw_3_Abbott[j] | log_beta_severity_Abbott[1], tau_Abbott[1]),
			normal_lpdf(log_intercept_raw_3_Abbott[j] | log_beta_severity_Abbott[2], tau_Abbott[2]));
		
		target += log_mix(theta_3_Roche,
			normal_lpdf(log_intercept_raw_3_Roche[j] | log_beta_severity_Roche[1], tau_Roche[1]),
			normal_lpdf(log_intercept_raw_3_Roche[j] | log_beta_severity_Roche[2], tau_Roche[2]));
		
	}
	
	// Likelihood part 2: data
	for (n in 1:N) {
		
		if(severity_by_obs[n]==1) {
			
			vector[2] mu_1;
			mu_1[1] = log_intercept_raw_1_Abbott[patient_ID_index_by_severity[n]] - (lambda_Abbott * time_since[n]);
			mu_1[2] = log_intercept_raw_1_Roche[patient_ID_index_by_severity[n]] - (lambda_Roche * time_since[n]);
			
			target += multi_normal_lpdf(log_SC_value[n,] | mu_1, Sigma);
			
		}
		
		if(severity_by_obs[n]==2) {
			
			vector[2] mu_2;
			mu_2[1] = log_intercept_raw_2_Abbott[patient_ID_index_by_severity[n]] - (lambda_Abbott * time_since[n]);
			mu_2[2] = log_intercept_raw_2_Roche[patient_ID_index_by_severity[n]] - (lambda_Roche * time_since[n]);
			
			target += multi_normal_lpdf(log_SC_value[n,] | mu_2, Sigma);
			
		}
		
		if(severity_by_obs[n]==3) {
			
			vector[2] mu_3;
			mu_3[1] = log_intercept_raw_3_Abbott[patient_ID_index_by_severity[n]] - (lambda_Abbott * time_since[n]);
			mu_3[2] = log_intercept_raw_3_Roche[patient_ID_index_by_severity[n]] - (lambda_Roche * time_since[n]);
			
			target += multi_normal_lpdf(log_SC_value[n,] | mu_3, Sigma);
			
		}
		
	}
	
	// Priors
	target += normal_lpdf(lambda_Abbott | 0.0, 0.01);
	target += normal_lpdf(lambda_Roche | 0.0, 0.01);
	
	target += normal_lpdf(log_beta_severity_Abbott[1] | 0.0, 2.5);
	target += normal_lpdf(log_beta_severity_Abbott[2] | 0.0, 2.5);

	target += normal_lpdf(log_beta_severity_Roche[1] | 0.0, 2.5);
	target += normal_lpdf(log_beta_severity_Roche[2] | 0.0, 2.5);
	
	target += normal_lpdf(tau_backwards_Abbott[1] | 0.0, 1.0);
	target += normal_lpdf(tau_backwards_Abbott[2] | 0.0, 1.0);
	
	target += normal_lpdf(tau_backwards_Roche[1] | 0.0, 1.0);
	target += normal_lpdf(tau_backwards_Roche[2] | 0.0, 1.0);
	
	target += normal_lpdf(theta_1_Abbott | 0.5, 0.5);
	target += normal_lpdf(theta_1_Roche | 0.5, 0.5);
	
	target += normal_lpdf(theta_2_Abbott | 0.5, 0.5);
	target += normal_lpdf(theta_2_Roche | 0.5, 0.5);
	
	target += normal_lpdf(theta_3_Abbott | 0.5, 0.5);
	target += normal_lpdf(theta_3_Roche | 0.5, 0.5);
	
}
generated quantities {
	
	matrix<lower=0.0, upper=1.0>[12,(DAYS_FOR_SENS+1)] sensitivity_over_time;	// Sensitivity: 3 classes of severity x 4 possible outcomes (++, +-, -+, --)
	
	// Calculate time-varying sensitivity by symptoms
	for(k in 0:DAYS_FOR_SENS) {
		
		matrix[2,N_rep] which_group_1_low; 		// Abbott, Roche
		matrix[2,N_rep] which_group_2_low; 		// Abbott, Roche
		matrix[2,N_rep] which_group_3_low; 		// Abbott, Roche
		
		matrix[2,N_rep] intercept_new_1_tau; 	// Abbott, Roche
		matrix[2,N_rep] intercept_new_2_tau;	// Abbott, Roche
		
		matrix[2,N_rep] Y_tilde_1;				// Abbott, Roche
		matrix[2,N_rep] Y_tilde_2;				// Abbott, Roche
		matrix[2,N_rep] Y_tilde_3;				// Abbott, Roche
		
		matrix[2,N_rep] Y_tilde_binary_1;		// Abbott, Roche
		matrix[2,N_rep] Y_tilde_binary_2;		// Abbott, Roche
		matrix[2,N_rep] Y_tilde_binary_3;		// Abbott, Roche
		
		vector[2] mu_1;
		vector[2] mu_2;
		vector[2] mu_3;
		
		// Simulate new data
		for(x in 1:N_rep) {
			
			// Assign the observation for class 1 as a low or high responder
			which_group_1_low[1,x] = bernoulli_rng(theta_1_Abbott);
			which_group_1_low[2,x] = bernoulli_rng(theta_1_Roche);
			
			// Assign the observation for class 2 as a low or high responder
			which_group_2_low[1,x] = bernoulli_rng(theta_2_Abbott);
			which_group_2_low[2,x] = bernoulli_rng(theta_2_Roche);
			
			// Assign the observation for class 3 as a low or high responder
			which_group_3_low[1,x] = bernoulli_rng(theta_3_Abbott);
			which_group_3_low[2,x] = bernoulli_rng(theta_3_Roche);
			
			// Generate other random variables
			intercept_new_1_tau[1,x] = normal_rng(0.0, tau_Abbott[1]);
			intercept_new_1_tau[2,x] = normal_rng(0.0, tau_Roche[1]);
			
			intercept_new_2_tau[1,x] = normal_rng(0.0, tau_Abbott[2]);
			intercept_new_2_tau[2,x] = normal_rng(0.0, tau_Roche[2]);
			
			// Generate continuous Y
			mu_1[1] = which_group_1_low[1,x] == 1 ?
					log_beta_severity_Abbott[1] + intercept_new_1_tau[1,x] - (lambda_Abbott * (k * 1.0)) :
					log_beta_severity_Abbott[2] + intercept_new_2_tau[1,x] - (lambda_Abbott * (k * 1.0));
			
			mu_1[2] = which_group_1_low[2,x] == 1 ?
					log_beta_severity_Roche[1] + intercept_new_1_tau[2,x] - (lambda_Roche * (k * 1.0)) :
					log_beta_severity_Roche[2] + intercept_new_2_tau[2,x] - (lambda_Roche * (k * 1.0));
			
			Y_tilde_1[,x] = multi_normal_rng(mu_1, Sigma);

			mu_2[1] = which_group_2_low[1,x] == 1 ?
					log_beta_severity_Abbott[1] + intercept_new_1_tau[1,x] - (lambda_Abbott * (k * 1.0)) :
					log_beta_severity_Abbott[2] + intercept_new_2_tau[1,x] - (lambda_Abbott * (k * 1.0));
			
			mu_2[2] = which_group_2_low[2,x] == 1 ?
					log_beta_severity_Roche[1] + intercept_new_1_tau[2,x] - (lambda_Roche * (k * 1.0)) :
					log_beta_severity_Roche[2] + intercept_new_2_tau[2,x] - (lambda_Roche * (k * 1.0));
			
			Y_tilde_2[,x] = multi_normal_rng(mu_2, Sigma);
			
			mu_3[1] = which_group_3_low[1,x] == 1 ?
					log_beta_severity_Abbott[1] + intercept_new_1_tau[1,x] - (lambda_Abbott * (k * 1.0)) :
					log_beta_severity_Abbott[2] + intercept_new_2_tau[1,x] - (lambda_Abbott * (k * 1.0));
			
			mu_3[2] = which_group_3_low[2,x] == 1 ?
					log_beta_severity_Roche[1] + intercept_new_1_tau[2,x] - (lambda_Roche * (k * 1.0)) :
					log_beta_severity_Roche[2] + intercept_new_2_tau[2,x] - (lambda_Roche * (k * 1.0));
			
			Y_tilde_3[,x] = multi_normal_rng(mu_3, Sigma);
			
			// Transform Y to binary results based on cutoff
			Y_tilde_binary_1[1,x] = Y_tilde_1[1,x] >= CUTOFF[1] ? 1.0 : 0.0;
			Y_tilde_binary_2[1,x] = Y_tilde_2[1,x] >= CUTOFF[1] ? 1.0 : 0.0;
			Y_tilde_binary_3[1,x] = Y_tilde_3[1,x] >= CUTOFF[1] ? 1.0 : 0.0;
			
			Y_tilde_binary_1[2,x] = Y_tilde_1[2,x] >= CUTOFF[2] ? 1.0 : 0.0;
			Y_tilde_binary_2[2,x] = Y_tilde_2[2,x] >= CUTOFF[2] ? 1.0 : 0.0;
			Y_tilde_binary_3[2,x] = Y_tilde_3[2,x] >= CUTOFF[2] ? 1.0 : 0.0;
			
		}
		
		// Calculate sensitivity
		sensitivity_over_time[1,k+1] = dot_product(Y_tilde_binary_1[1,], Y_tilde_binary_1[2,]) / (1.0 * N_rep);
		sensitivity_over_time[2,k+1] = dot_product(Y_tilde_binary_1[1,], (1.0-Y_tilde_binary_1[2,])) / (1.0 * N_rep);
		sensitivity_over_time[3,k+1] = dot_product((1.0-Y_tilde_binary_1[1,]), Y_tilde_binary_1[2,]) / (1.0 * N_rep);
		sensitivity_over_time[4,k+1] = dot_product((1.0-Y_tilde_binary_1[1,]), (1.0-Y_tilde_binary_1[2,])) / (1.0 * N_rep);
		
		sensitivity_over_time[5,k+1] = dot_product(Y_tilde_binary_2[1,], Y_tilde_binary_2[2,]) / (1.0 * N_rep);
		sensitivity_over_time[6,k+1] = dot_product(Y_tilde_binary_2[1,], (1.0-Y_tilde_binary_2[2,])) / (1.0 * N_rep);
		sensitivity_over_time[7,k+1] = dot_product((1.0-Y_tilde_binary_2[1,]), Y_tilde_binary_2[2,]) / (1.0 * N_rep);
		sensitivity_over_time[8,k+1] = dot_product((1.0-Y_tilde_binary_2[1,]), (1.0-Y_tilde_binary_2[2,])) / (1.0 * N_rep);
		
		sensitivity_over_time[9,k+1] = dot_product(Y_tilde_binary_3[1,], Y_tilde_binary_3[2,]) / (1.0 * N_rep);
		sensitivity_over_time[10,k+1] = dot_product(Y_tilde_binary_3[1,], (1.0-Y_tilde_binary_3[2,])) / (1.0 * N_rep);
		sensitivity_over_time[11,k+1] = dot_product((1.0-Y_tilde_binary_3[1,]), Y_tilde_binary_3[2,]) / (1.0 * N_rep);
		sensitivity_over_time[12,k+1] = dot_product((1.0-Y_tilde_binary_3[1,]), (1.0-Y_tilde_binary_3[2,])) / (1.0 * N_rep);
		
	}
	
}
