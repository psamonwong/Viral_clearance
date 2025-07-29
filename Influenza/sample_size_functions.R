# function that simulates one patient given parameter settings
# this assumes log-linear decline in viral loads

f_sim = function(t_design, # design sampling times for PCR swabs
                 my_intercept, # intercept
                 my_slope, # slope
                 sigma_vl, # standard deviation for t-distribution error model
                 t_dof, # degrees of freedom for the t-distribution error model
                 LOD # lower limit of quantification/detection
){

  true_log_vl = my_intercept + my_slope*t_design

  # Add error using a t-distribution
  log_vl_sim =
    LaplacesDemon::rst(n = length(t_design),
                       mu = true_log_vl,
                       sigma = sigma_vl,
                       nu = t_dof)

  if(any(is.na(log_vl_sim))) {
    print('generated NA values!!!')
  }

  # Values below LOD set to LOD
  log_vl_sim = ifelse(log_vl_sim<LOD, LOD, log_vl_sim)

  out_sim = data.frame(log_vl_sim=log_vl_sim,
                       true_log_vl=true_log_vl)
  return(out_sim)

}

# wrapper function
sim_individuals = function(thetas, # posterior distribution: a stan object
                           t_design, # design points for the swabs
                           N, # total number of patients
                           trt_effects, # effects for each arm
                           Trt_arm, # vector of length N with treatment arms for each patient
                           LOD = 1, # lower limit of detection for censoring
                           f_sim
){

  if(!(N == length(Trt_arm))) stop('Trt_arm has to be of length=N')
  if(!(length(unique(Trt_arm)) == length(trt_effects))) stop('trt_effects has to be of length number of unique Trt_Arm')
  
  K=length(thetas$t_dof) # number of posterior samples in model fit
  post_i = sample(x = K, size = 1) # choose a random posterior draw on which to base simulation

  # make the variance/covariance matrix for the random effects
  L = diag(x = unlist(thetas$sigmasq_u[post_i,]),nrow = 2,ncol = 2)%*%thetas$L_Omega[post_i,,]
  Epsilon = L%*%t(L)

  # generate random effects matrix for the N patients
  thetas_rand = mvtnorm::rmvnorm(n = N, sigma = Epsilon)

  # make sim data matrix
  Log_VL = array(dim = c(length(t_design)*N, 7))
  Log_VL = as.data.frame(Log_VL)
  colnames(Log_VL) = c('log10_viral_load',
                       'log10_cens_vl',
                       'True_Log_VL',
                       'Time',
                       'ID',
                       'Trt_effect',
                       'Trt_arm')
  Log_VL$ID = as.vector(sapply(1:N, rep, length(t_design)))
  Log_VL$log10_cens_vl=LOD
  # make sure the t_design points are sorted in increasing order
  t_design = sort(t_design)
  for(n in 1:N){

    ind_patient = which(Log_VL$ID==n)
    xs = f_sim(t_design = t_design,
               my_intercept = thetas$alpha_0[post_i]+thetas_rand[n,1],
               my_slope = thetas$beta_0[post_i] * exp(thetas_rand[n,2] + log(trt_effects[Trt_arm[n]])),
               sigma_vl = thetas$sigma_logvl[post_i],
               t_dof = thetas$t_dof[post_i],
               LOD=LOD)



    Log_VL[ind_patient, 'log10_viral_load'] = xs$log_vl_sim
    Log_VL[ind_patient, 'True_Log_VL'] = xs$true_log_vl

    Log_VL[ind_patient, 'Time'] = t_design
    Log_VL[ind_patient, 'Trt_effect'] = trt_effects[Trt_arm[n]]
    Log_VL[ind_patient, 'Trt_arm'] = Trt_arm[n]
    
  }
  return(Log_VL)
}


make_stan_inputs = function(input_data_fit,
                            trt_frmla,
                            Dmax
){

  ## check censored values come last
  if(!all(diff(input_data_fit$log10_viral_load == input_data_fit$log10_cens_vl)>=0)) stop('problem with ordering')

  ind_dup = !duplicated(input_data_fit$ID)

  ind_cens = !input_data_fit$log10_viral_load>
    input_data_fit$log10_cens_vl
  ID_map = data.frame(ID_key = input_data_fit$ID,
                      ID_stan = as.numeric(as.factor(input_data_fit$ID)))
  analysis_data_stan = list(Ntot = nrow(input_data_fit),
                            N_obs = sum(!ind_cens),
                            n_id = max(ID_map$ID_stan),
                            id = ID_map$ID_stan,
                            ind_start = which(!duplicated(ID_map$ID_stan)),
                            obs_day = input_data_fit$Time,
                            log_10_vl = input_data_fit$log10_viral_load,
                            log10_cens_vl = input_data_fit$log10_cens_vl,
                            Time_max = Dmax)
  ID_map = ID_map[!duplicated(ID_map$ID_key), ]

  all(analysis_data_stan$log_10_vl[1:analysis_data_stan$N_obs]>
        analysis_data_stan$log10_cens_vl[1:analysis_data_stan$N_obs]) &
    all(analysis_data_stan$log_10_vl[(1+analysis_data_stan$N_obs):analysis_data_stan$Ntot] ==
          analysis_data_stan$log10_cens_vl[(1+analysis_data_stan$N_obs):analysis_data_stan$Ntot])

  Trt_matrix = model.matrix(trt_frmla, data = input_data_fit)
  Trt_matrix[,1]=0 # first column is dummy

  analysis_data_stan$trt_mat = Trt_matrix
  analysis_data_stan$K_trt = ncol(Trt_matrix)

  return(analysis_data_stan)
}

# Checks a function for use of global variables
# Returns TRUE if ok, FALSE if globals were found.
checkStrict <- function(f, silent=FALSE) {
  vars <- codetools::findGlobals(f)
  found <- !vapply(vars, exists, logical(1), envir=as.environment(2))
  if (!silent && any(found)) {
    warning("global variables used: ", paste(names(found)[found], collapse=', '))
    return(invisible(FALSE))
  }

  !any(found)
}



# wrapper function
sim_individuals_biexp = function(thetas, # posterior distribution: a stan object
                                 t_design, # design points for the swabs
                                 N, # total number of patients
                                 trt_effects, # effects for each arm
                                 Trt_arm, # vector of length N with treatment arms for each patient
                                 LOD = 1, # lower limit of detection for censoring
                                 f_sim_biexp
){
  
  if(!(N == length(Trt_arm))) stop('Trt_arm has to be of length=N')
  if(!(length(unique(Trt_arm)) == length(trt_effects))) stop('trt_effects has to be of length number of unique Trt_Arm')
  
  K=length(thetas$t_dof) # number of posterior samples in model fit
  post_i = sample(x = K, size = 1) # choose a random posterior draw on which to base simulation
  
  # make the variance/covariance matrix for the random effects
  L = diag(x = unlist(thetas$sigmasq_u[post_i,]),nrow = 4,ncol = 4)%*%thetas$L_Omega[post_i,,]
  Epsilon = L%*%t(L)
  
  # generate random effects matrix for the N patients
  thetas_rand = mvtnorm::rmvnorm(n = N, sigma = Epsilon)
  
  # make sim data matrix
  Log_VL = array(dim = c(length(t_design)*N, 7))
  Log_VL = as.data.frame(Log_VL)
  colnames(Log_VL) = c('log10_viral_load',
                       'log10_cens_vl',
                       'True_Log_VL',
                       'Time',
                       'ID',
                       'Trt_effect',
                       'Trt_arm')
  Log_VL$ID = as.vector(sapply(1:N, rep, length(t_design)))
  Log_VL$log10_cens_vl=LOD
  # make sure the t_design points are sorted in increasing order
  t_design = sort(t_design)
  for(n in 1:N){
    
    ind_patient = which(Log_VL$ID==n)
    xs = f_sim_biexp(t_design = t_design,
                     my_intercept_1 = thetas$intercept[post_i,1]+thetas_rand[n,3],
                     my_intercept_2 = thetas$intercept[post_i,2]+thetas_rand[n,1],
                     my_slope_1 = thetas$coef[post_i,1] * exp(thetas_rand[n,4]),
                     my_slope_2 = thetas$coef[post_i,2] * exp(thetas_rand[n,2] + log(trt_effects[Trt_arm[n]])),
                     #my_intercept = thetas$alpha_0[post_i]+thetas_rand[n,1],
                     #my_slope = thetas$beta_0[post_i] * exp(thetas_rand[n,2] + log(trt_effects[Trt_arm[n]])),
                     sigma_vl = thetas$sigma_logvl[post_i],
                     t_dof = thetas$t_dof[post_i],
                     LOD=LOD)
    
    
    
    Log_VL[ind_patient, 'log10_viral_load'] = xs$log_vl_sim
    Log_VL[ind_patient, 'True_Log_VL'] = xs$true_log_vl
    
    Log_VL[ind_patient, 'Time'] = t_design
    Log_VL[ind_patient, 'Trt_effect'] = trt_effects[Trt_arm[n]]
    Log_VL[ind_patient, 'Trt_arm'] = Trt_arm[n]
    
  }
  return(Log_VL)
}






library(matrixStats)


# function that simulates one patient given parameter settings
# this assumes bi-exponential decline in viral loads

f_sim_biexp = function(t_design, # design sampling times for PCR swabs
                       my_intercept_1, # intercept1
                       my_intercept_2, # intercept2
                       my_slope_1, # slope1
                       my_slope_2, # slope2
                       sigma_vl, # standard deviation for t-distribution error model
                       t_dof, # degrees of freedom for the t-distribution error model
                       LOD # lower limit of quantification/detection
){
  
  #true_log_vl = my_intercept + my_slope*t_design
  
  # Vectorized computation
  term1 <- my_intercept_2 - my_slope_2 * t_design
  term2 <- my_intercept_1 - my_slope_1 * t_design
  
  # Combine terms into matrix for logSumExp
  terms_matrix <- rbind(term1, term2)
  
  # Apply logSumExp column-wise
  true_log_vl <- apply(terms_matrix, 2, matrixStats::logSumExp)
  
  # Add error using a t-distribution
  log_vl_sim =
    LaplacesDemon::rst(n = length(t_design),
                       mu = true_log_vl,
                       sigma = sigma_vl,
                       nu = t_dof)
  
  if(any(is.na(log_vl_sim))) {
    print('generated NA values!!!')
  }
  
  # Values below LOD set to LOD
  log_vl_sim = ifelse(log_vl_sim<LOD, LOD, log_vl_sim)
  
  out_sim = data.frame(log_vl_sim=log_vl_sim,
                       true_log_vl=true_log_vl)
  return(out_sim)
  
}


# wrapper function
sim_individuals_growth_decay= function(thetas, # posterior distribution: a stan object
                                       t_design, # design points for the swabs
                                       N, # total number of patients
                                       trt_effects, # effects for each arm
                                       Trt_arm, # vector of length N with treatment arms for each patient
                                       LOD = 1, # lower limit of detection for censoring
                                       f_sim_growth_decay
){
  
  if(!(N == length(Trt_arm))) stop('Trt_arm has to be of length=N')
  if(!(length(unique(Trt_arm)) == length(trt_effects))) stop('trt_effects has to be of length number of unique Trt_Arm')
  
  K=length(thetas$t_dof) # number of posterior samples in model fit
  post_i = sample(x = K, size = 1) # choose a random posterior draw on which to base simulation
  
  # make the variance/covariance matrix for the random effects
  L = diag(x = unlist(thetas$sigmasq_u[post_i,]),nrow = 3,ncol = 3)%*%thetas$L_Omega[post_i,,]
  Epsilon = L%*%t(L)
  
  # generate random effects matrix for the N patients
  thetas_rand = mvtnorm::rmvnorm(n = N, sigma = Epsilon)
  
  # make sim data matrix
  Log_VL = array(dim = c(length(t_design)*N, 7))
  Log_VL = as.data.frame(Log_VL)
  colnames(Log_VL) = c('log10_viral_load',
                       'log10_cens_vl',
                       'True_Log_VL',
                       'Time',
                       'ID',
                       'Trt_effect',
                       'Trt_arm')
  Log_VL$ID = as.vector(sapply(1:N, rep, length(t_design)))
  Log_VL$log10_cens_vl=LOD
  # make sure the t_design points are sorted in increasing order
  t_design = sort(t_design)
  for(n in 1:N){
    
    ind_patient = which(Log_VL$ID==n)
    xs = f_sim_growth_decay(t_design = t_design,
                            my_intercept = thetas$alpha_0[post_i]+thetas_rand[n,1],
                            my_a = thetas$beta_0[post_i,1],
                            my_b = thetas$beta_0[post_i,2] * exp(thetas_rand[n,2] + log(trt_effects[Trt_arm[n]])),
                            my_tmax = thetas$tmax_pop[post_i] + thetas_rand[n,3],
                            sigma_vl = thetas$sigma_logvl[post_i],
                            t_dof = thetas$t_dof[post_i],
                            LOD=LOD)
    
    Log_VL[ind_patient, 'log10_viral_load'] = xs$log_vl_sim
    Log_VL[ind_patient, 'True_Log_VL'] = xs$true_log_vl
    
    Log_VL[ind_patient, 'Time'] = t_design
    Log_VL[ind_patient, 'Trt_effect'] = trt_effects[Trt_arm[n]]
    Log_VL[ind_patient, 'Trt_arm'] = Trt_arm[n]
    
  }
  return(Log_VL)
}


f_sim_growth_decay = function(t_design, # design sampling times for PCR swabs
                              my_intercept,
                              my_a,
                              my_b,
                              my_tmax,
                              sigma_vl, # standard deviation for t-distribution error model
                              t_dof, # degrees of freedom for the t-distribution error model
                              LOD # lower limit of quantification/detection
){
  
  # Vectorized calculation
  term1 <- log(my_b) - my_a * (t_design - my_tmax)
  term2 <- log(my_a) + my_b * (t_design - my_tmax)
  
  # Combine terms into matrix for logSumExp
  terms_matrix <- rbind(term1, term2)
  
  # Apply logSumExp column-wise and compute final prediction
  true_log_vl <- my_intercept + log(my_a + my_b) - apply(terms_matrix, 2, matrixStats::logSumExp)
  
  
  # Add error using a t-distribution
  log_vl_sim =
    LaplacesDemon::rst(n = length(t_design),
                       mu = true_log_vl,
                       sigma = sigma_vl,
                       nu = t_dof)
  
  if(any(is.na(log_vl_sim))) {
    print('generated NA values!!!')
  }
  
  # Values below LOD set to LOD
  log_vl_sim = ifelse(log_vl_sim<LOD, LOD, log_vl_sim)
  
  out_sim = data.frame(log_vl_sim=log_vl_sim,
                       true_log_vl=true_log_vl)
  return(out_sim)
  
}


checkStrict(make_stan_inputs)
checkStrict(sim_individuals)
checkStrict(f_sim)

