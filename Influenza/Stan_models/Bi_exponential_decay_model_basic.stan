//   Copyright 2024 James Watson and Phrutsamon Wongnak, Mahidol Oxford Tropical Medicine Research Unit
//
//   Licensed under a CC-BY License
/**
Model 0 used for the simulation studies with no treatment effects and RNaseP
- Bi-exponential decay
**/


data {
  int<lower=0> Ntot;                           // Number of PCR data points
  int<lower=0,upper=Ntot> N_obs;               // Number of PCR data points with viral load less than LOD
  int<lower=0> n_id;                           // Number of individuals
  int<lower=1,upper=n_id> id[Ntot];            // Patient identifier for each PCR sample
  int<lower=1,upper=Ntot> ind_start[n_id];     // Starting index for each patient: used in generated quantities block
  real<lower=0> Time_max;                      // Max follow-up time
  real<lower=0,upper=Time_max> obs_day[Ntot];                 // Vector marking the time for each data point
  
  real log_10_vl[Ntot];                        // log base 10 viral load in copies per mL
  real log10_cens_vl[Ntot];                    // censoring value for censored observation

  // priors
  real A0_prior;
  real B0_prior;
  real coef_1_prior;
  real coef_2_prior;
  
  real<lower=0> prior_intercept_sd;
  real<lower=0> prior_coef_sd;
}

transformed data {
  vector[4] zeros;
  for(i in 1:4) zeros[i] = 0;
}

parameters {
  // hyperparameters
  cholesky_factor_corr[4] L_Omega;      // correlation matrix
  vector<lower=0>[4] sigmasq_u;         // variance of random effects
  
  real<lower=0> sigma_vl;
  real<lower=0> t_dof;                  // student-t degrees of freedom

  ordered[2] intercept;                 // intercept for second and first components respectively
  ordered[2] coef;                      // slope for second and first components respectively
  vector[4] theta_rand[n_id];           // individual random effects vector
  
}

transformed parameters {
  real pred_log10_vl[Ntot];
  for(i in 1:Ntot){
    pred_log10_vl[i] = log_sum_exp(intercept[2]+theta_rand[id[i]][1]-
    (coef[2]*exp(theta_rand[id[i]][2]))*obs_day[i],
    intercept[1]+theta_rand[id[i]][3]-
    (coef[1]*exp(theta_rand[id[i]][4]))*obs_day[i]); 
  }
}

model {
  // random effects
  sigmasq_u[1] ~ exponential(1);
  sigmasq_u[2] ~ exponential(1);
  sigmasq_u[3] ~ exponential(1);
  sigmasq_u[4] ~ exponential(1);
  
  // measurement error
  sigma_vl ~ exponential(1);
  t_dof ~ exponential(1);

  // covariance matrix - random effects
  L_Omega ~ lkj_corr_cholesky(1);
  for (i in 1:n_id) theta_rand[i] ~ multi_normal_cholesky(zeros, diag_pre_multiply(sigmasq_u, L_Omega));  
  
  intercept[2] ~ normal(A0_prior,prior_intercept_sd);
  intercept[1] ~ normal(B0_prior,prior_intercept_sd);
  coef[2] ~ normal(coef_1_prior, prior_coef_sd);
  coef[1] ~ normal(coef_2_prior, prior_coef_sd);
  
  //***** Likelihood *****
  // Non censored observations
  log_10_vl[1:N_obs] ~ student_t(t_dof, pred_log10_vl[1:N_obs], sigma_vl);
  
  // Censored observations
  for(i in (N_obs+1):Ntot){
    target += student_t_lcdf(log10_cens_vl | t_dof, pred_log10_vl[i], sigma_vl);
  }
}

generated quantities {
  real preds[Ntot]; // For plotting
  vector[n_id] slope_1;
  vector[n_id] slope_2;

  
  for(i in 1:N_obs){
    preds[i] = pred_log10_vl[i];
  }
  for(i in (N_obs+1):Ntot){
    preds[i] = pred_log10_vl[i];
  }
  
  for(i in 1:n_id){
    int j = ind_start[i];
    slope_1[i] = (coef[1]*exp(theta_rand[id[j]][4]));
    slope_2[i] = (coef[2]*exp(theta_rand[id[j]][2]));
  }
}


