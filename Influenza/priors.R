get_priors = function(){
  prior_params = list(alpha_0_prior_mean = 5,
                      alpha_0_prior_sd = 2,
                      beta_0_prior_mean = -.5,
                      beta_0_prior_sd = 1,
                      trt_effect_sd = 1,
                      sigma_logvl_mean = 1,
                      sigma_logvl_sd = 1,
                      gamma_rnasep_sd = 1,
                      slope_coefs_sd = .5,
                      intercept_coefs_sd = .5,
                      t_max_pop_mean = -3,
                      t_max_pop_sd = 3,
                      A0_prior = 5,
                      B0_prior = 3,
                      coef_1_prior = 0.25,
                      coef_2_prior = 1,
                      prior_intercept_sd = 2,
                      prior_coef_sd = 1)
  
  prior_params_sensitive = prior_params
  for(i in grep('sd',x = names(prior_params))){
    prior_params_sensitive[[i]] = prior_params_sensitive[[i]]*5
  }
  
  ## we get convergence issues with the covariate models if you have very wide priors on the slope coefs
  prior_params_sensitive$slope_coefs_sd=1
  
  all_priors = list(WIP=prior_params, NIP=prior_params_sensitive)
  
  return(all_priors)
}
all_priors=get_priors()
