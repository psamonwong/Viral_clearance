library(dplyr)
library(tidyr)
############################################################################

############################################################################
source('priors.R')
source('functions.R')
############################################################################
mITT <- function(platcov_dat){
  platcov_dat = platcov_dat %>% group_by(ID) %>%
    mutate(
      mITT = any(Per_protocol_sample==1 & Timepoint_ID>=3) & !all(CT_NS==40))
  return(as.data.frame(platcov_dat))
}

############################################################################
# Loading data
data = read.csv('Analysis_Data/no_study_drug_analysis.csv')

data %>%
  mutate(count = 1) %>%
  group_by(ID, Timepoint_ID) %>%
  


data <- mITT(recent_paxlovid_data)




data_list <- list(#remdesivir_data, 
                  #paxlovid_data, 
                  recent_paxlovid_data#, 
                  #ineffective_data
                  )

############################################################################
interventions_all = c("Nirmatrelvir + Ritonavir")
ref_arms_all = c("No study drug")
data_ID <- c(1)
pairs <- paste(interventions_all,ref_arms_all, data_ID, sep = "_")

Dmax_all <- 5 #c(2:7, 14)

#bootstrap_rep = 20

#pair_ID <- 6:8 # For ineffective drugs (should be 1:5 for success drugs)

model_settings <-  expand.grid(mod = 'Linear_model_basic.stan',
                               prior = 1,
                               cov_matrices = 1,
                               Dmax = Dmax_all,
                               pairs = pairs)#[pair_ID],
                               #boot_rep = 1:bootstrap_rep)

model_settings <- model_settings %>%
  separate(pairs, c("intervention", "ref_arm", "data_ID"), "_")

model_settings$Niter = 4000
model_settings$Nwarmup = 2000
model_settings$Nthin = 8
model_settings$Nchain = 4

writeLines(sprintf('Numbers of posterior samples in total is %s',
                   unique((model_settings$Niter-model_settings$Nwarmup)*model_settings$Nchain/model_settings$Nthin)))
model_settings$data_ID <- as.numeric(model_settings$data_ID)

############################################################################
model_setup_f_name = "model_settings.RData"
save(model_settings, 
     data_list,
     all_priors,
     file = model_setup_f_name)
