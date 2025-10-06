# 1. Load required libraries
# Load core tidyverse for data manipulation and plotting, and here for file path management
library(tidyverse)
library(here)

# 2. Define constants
# Set quantiles used for summarizing simulation results
QUANTILES <- c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
rerun <- FALSE

# Summarize simulation results: quantiles, mean, and sd
summarize_simulation <- function(sim_df, quantiles) {
  df_summary <- sim_df %>%
    reframe(across(everything(), ~ quantile(.x, quantiles))) %>%
    t() %>%
    as_tibble() 
  df_mean <- sim_df %>%
    reframe(across(everything(), ~ mean(.x))) %>%
    t() %>%
    as_tibble() 
  df_sd <- sim_df %>%
    reframe(across(everything(), ~ sd(.x))) %>%
    t() %>%
    as_tibble() 
  list(summary = df_summary, mean = df_mean, sd = df_sd)
}


# 4. Load data
# Load simulation settings and simulation results
load_simulation <- function(rdata, simdir, neg_control,  rerun, output_dir){
  load(here('Rout', rdata))
  sim_dir <- here(simdir)
  sim_settings <- sim_settings %>% mutate(i = row_number())
  if(neg_control){sim_settings <- sim_settings %>% filter(trt_control == 1)} else {sim_settings <- sim_settings %>% filter(trt_control != 1)}
  
  sim_list <- list.files(sim_dir, full.names = TRUE)
  sim_ids <- str_extract(sim_list, "out[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric()
  
  if (rerun) {
    # If rerun is TRUE, re-aggregate simulation results from individual files
    sim_tbl <- tibble(file = sim_list, id = sim_ids)
    tbl_results <- sim_settings %>%
      left_join(sim_tbl, by = c("i" = "id")) %>%
      mutate(res = map(file, ~ read_csv(.x, show_col_types = FALSE)$trt_effect))
    simulation_results <- tbl_results %>%
      pull(res) %>%
      set_names(paste0("sim_", seq_along(.))) %>%
      as_tibble() 
    write_csv(simulation_results, here("Outputs", output_dir))
  } else {
    # Otherwise, load the pre-aggregated simulation results
    simulation_results <- read_csv(here("Outputs", output_dir))
  }
  
  summary_stats <- summarize_simulation(simulation_results, QUANTILES)
  summary_quantiles <- summary_stats$summary
  mean_stats <- summary_stats$mean
  sd_stats <- summary_stats$sd
  colnames(summary_quantiles) <- c("low95", "low90", "q1", "q2", 'q3', "up90", "up95")
  
  all_summary_stats <- cbind(sim_settings, summary_quantiles, mean = mean_stats$V1, sd = sd_stats$V1)
  all_summary_stats <- all_summary_stats %>% mutate(error =  mean - log(trt_effect_comp))
  all_summary_stats
}


# Calculate power table for each scenario
calculate_power_table <- function(df_summary_all) {
  df_summary_all %>%
    group_by(day_plans, N, N_swabs_per_day, trt_effect_comp, k_slope, k_sigma_logvl, k_sigmasq_u_1,
             k_sigmasq_u_2, trt_control) %>%
    mutate(rejectH0 = low95 > 0) %>%
    summarise(power = mean(rejectH0), n = n(), .groups = 'drop') #%>%
    # relabel_factors() %>%
    # assign_highlight_type()
}

# 6. Calculate power and prepare data for modeling
# Calculate power for each scenario and prepare data for GAM fitting
power_summary <- calculate_power_table(all_summary_stats)

ggplot(power_summary %>% filter(N_swabs_per_day == 2,
                                #day_plans == "0,1,2,3,4,5",
                                trt_control == 1,
                                k_slope == 1,
                                k_sigma_logvl == 1,
                                k_sigmasq_u_1 == 1,
                                k_sigmasq_u_2 == 1),
       aes(x = N, y = power, group = day_plans, col = day_plans)) +
  geom_point()  +
  geom_line() +
  scale_x_log10() +
  theme_bw() +
  facet_wrap(.~as.factor(trt_effect_comp),ncol = 3)

prepare_gam_data <- function(df_summary_all) {
  df_summary_all %>%
    group_by(day_plans, N, N_swabs_per_day, trt_effect_comp, 
             k_slope, k_sigma_logvl, k_sigmasq_u_1, k_sigmasq_u_2
             ) %>%
    mutate(rejectH0 = low95 > 0) 
}


# Fit a GAM model and predict power curves
fit_gam_and_predict <- function(dat) {
  dat$trt_effect_comp <- as.factor(dat$trt_effect_comp)

  
  
  fit <- mgcv::gam(rejectH0 ~ log(N)*day_plans*trt_effect_comp, data = dat, family = "binomial")
  new_dat <- expand.grid(N = seq(10,2000,1),
                         day_plans = levels(dat$day_plans),
                         trt_effect_comp =  levels(as.factor(dat$trt_effect_comp)),
                         N_swabs_per_day =  levels(as.factor(dat$N_swabs_per_day)))
 # new_dat <- assign_highlight_type(new_dat)
  new_dat$preds = predict(fit, newdata = new_dat, type = 'response') %>% as.numeric()
  list(fit = fit, new_dat = new_dat)
}





all_summary_stats <- load_simulation(rdata = 'sim_settings_Jun25.RData',
                simdir = 'sims_out_June2025',
                neg_control = TRUE,
                rerun = TRUE,
                output_dir = "sim_posterior_Jun25_negative_control.csv") 
gam_data <- prepare_gam_data(all_summary_stats)


all_summary_stats %>% glimpse()

dat <- gam_data
dat$trt_effect_comp <- as.factor(dat$trt_effect_comp)

dat_for_fit <- dat %>% distinct(day_plans, N, N_swabs_per_day, trt_effect_comp, 
                                k_slope, k_sigma_logvl, k_sigmasq_u_1, k_sigmasq_u_2, .keep_all = T)

fit <- mgcv::gam(rejectH0 ~ log(N)*day_plans*trt_effect_comp*N_swabs_per_day*k_slope*k_sigma_logvl*k_sigmasq_u_1*k_sigmasq_u_2, 
                 data = dat, family = "binomial")










gam_fit_results <- fit_gam_and_predict(gam_data)
predicted_power <- gam_fit_results$new_dat
error_summary <- summarize_error_table(all_summary_stats)


ggplot(new_dat, aes(x = N, y = preds, col = day_plans)) +
  geom_line() +
  facet_wrap(.~as.factor(trt_effect_comp), ncol = 3) +
  #scale_x_log10() +
  theme_bw(base_size = 15)