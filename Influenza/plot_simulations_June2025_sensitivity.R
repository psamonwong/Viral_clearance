# ==============================================================================
# SENSITIVITY ANALYSIS FOR VIRAL CLEARANCE POWER CALCULATIONS
# ==============================================================================
# This script compares power curves across different data generating processes:
# 1) Exponential decay model (baseline)
# 2) Bi-exponential decay model 
# 3) Growth-and-decay model
# ==============================================================================

# 1. SETUP & DEPENDENCIES =====================================================
# Load required libraries
library(tidyverse)    # Data manipulation and plotting
library(here)         # File path management
library(mgcv)         # GAM modeling

# 2. CONFIGURATION ============================================================
# Set quantiles used for summarizing simulation results
QUANTILES <- c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
rerun <- FALSE
power_threshold <- 0.8

# Define plot styling constants
TREATMENT_EFFECT_LABELS <- c("Treatment effect = 0%", "Treatment effect = 20%", 
                            "Treatment effect = 40%", "Treatment effect = 60%", 
                            "Treatment effect = 80%", "Treatment effect = 100%")
COLOR_PALETTE <- rev(c("#1230AE", "#6C48C5", "#C68FE6", "#D8A25E", "#A04747", "black"))
N_BREAKS <- c(10, 20, 40, 80, 160, 320, 640)

# 3. HELPER FUNCTIONS ==========================================================

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

# Load simulation settings and simulation results
load_simulation <- function(rdata, simdir, neg_control, rerun, output_dir){
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

# Prepare data for GAM modeling with model-specific grouping
prepare_gam_data <- function(df_summary_all, baseline) {
  if(baseline){ 
    df_summary_all %>%
      group_by(day_plans, N, N_swabs_per_day, trt_effect_comp, 
               k_slope, k_sigma_logvl, k_sigmasq_u_1, k_sigmasq_u_2) %>%
      mutate(rejectH0 = low95 > 0)
  } else {
    df_summary_all %>%
      group_by(day_plans, N, N_swabs_per_day, trt_effect_comp) %>%
      mutate(rejectH0 = low95 > 0)  
  }
}

# 4. DATA LOADING ==============================================================

# Load exponential decay model (baseline) results
all_summary_stats_linear <- load_simulation(
  rdata = 'sim_settings_Jun25.RData',
  simdir = 'sims_out_June2025',
  neg_control = TRUE,
  rerun = F,
  output_dir = "sim_posterior_Jun25_negative_control.csv"
) 

# Filter to baseline conditions and prepare data
gam_data_linear <- prepare_gam_data(all_summary_stats_linear, T)
gam_data_linear <- gam_data_linear %>% 
  filter(N_swabs_per_day == 2,
         day_plans == "0,1,2,3,4,5",
         k_slope == 1,
         k_sigma_logvl == 1,
         k_sigmasq_u_1 == 1,
         k_sigmasq_u_2 == 1) %>%
  select(-k_slope, -k_sigma_logvl, -k_sigmasq_u_1, -k_sigmasq_u_2)
gam_data_linear$mod <- "Exponential decay model"

# Load bi-exponential decay model results
all_summary_stats_biexp <- load_simulation(
  rdata = 'sim_settings_Jun25_biexponential.RData',
  simdir = 'sims_out_June2025_biexp',
  neg_control = TRUE,
  rerun = F,
  output_dir = "sim_posterior_Jun25_biexponential.csv"
) 
gam_data_biexp <- prepare_gam_data(all_summary_stats_biexp, F)
gam_data_biexp$mod <- "Bi-exponential decay model"

# Load growth-and-decay model results
all_summary_stats_growth_decay <- load_simulation(
  rdata = 'sim_settings_Jun25_growth_decay.RData',
  simdir = 'sims_out_June2025_growth_decay',
  neg_control = TRUE,
  rerun = F,
  output_dir = "sim_posterior_Jun25_growth_decay.csv"
) 
gam_data_growth_decay <- prepare_gam_data(all_summary_stats_growth_decay, F)
gam_data_growth_decay$mod <- "Growth-and-decay model"

# Combine all model results
gam_data <- rbind(gam_data_linear, gam_data_biexp, gam_data_growth_decay)

# 5. SENSITIVITY ANALYSIS PLOT =================================================

# Prepare data for plotting
dat_for_fit <- gam_data

# Format treatment effect and model labels
dat_for_fit$trt_effect_comp <- as.factor(dat_for_fit$trt_effect_comp)
levels(dat_for_fit$trt_effect_comp) <- TREATMENT_EFFECT_LABELS

dat_for_fit$mod <- factor(dat_for_fit$mod,
                         levels = c("Exponential decay model", 
                                   "Bi-exponential decay model", 
                                   "Growth-and-decay model"))

# Calculate empirical power at each sample size by model
summary_p <- dat_for_fit %>%
  group_by(N, trt_effect_comp, mod) %>%
  summarise(n = n(), p_rejectH0 = sum(rejectH0) / n, .groups = "drop")

# Fit GAM model for smooth power curves across models
fit <- mgcv::gam(rejectH0 ~ log10(N)*trt_effect_comp*mod, 
                 data = dat_for_fit, family = "binomial")

# Generate predictions for smooth curves
new_dat <- expand.grid(
  N = seq(1, 2000, 1),
  trt_effect_comp = levels(as.factor(dat_for_fit$trt_effect_comp)),
  mod = levels(as.factor(dat_for_fit$mod))
)
new_dat$preds <- predict(fit, newdata = new_dat, type = 'response') %>% as.numeric()

# Create sensitivity analysis plot
sensitivity_plot <- ggplot(mapping = aes(col = trt_effect_comp)) +
  geom_line(data = new_dat, mapping = aes(x = N, y = preds),
            linewidth = 1, alpha = 0.75, linetype = '31') +
  geom_point(data = summary_p, mapping = aes(x = N, y = p_rejectH0),
             size = 3, alpha = 0.75) +
  scale_x_log10(breaks = N_BREAKS, limits = c(10, 640)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 15) +
  facet_wrap(.~mod, ncol = 3) +
  scale_color_manual(values = COLOR_PALETTE, name = "") +
  geom_hline(yintercept = c(0.8, 0.9), col = "gray20", linetype = "22", linewidth = 0.75, alpha = 0.8) +
  geom_hline(yintercept = c(0, 1), col = "grey20", linetype = "22", 
             linewidth = 0.75, alpha = 0.2) +
  xlab("Number of patients per arm") +
  ylab("Pr[Treatment effect ≤ 0] < 0.025") +
  ggtitle("Sensitivity analysis for data generating processes") 

# Display plot
print(sensitivity_plot)

# 6. SAMPLE SIZE CALCULATIONS ==================================================

# Calculate sample sizes needed for 80% power by model
sample_sizes_by_model <- new_dat %>%
  mutate(diff_pred = abs(preds - power_threshold)) %>%
  group_by(trt_effect_comp, mod) %>%
  filter(diff_pred == min(diff_pred) & diff_pred < 0.01) %>%
  ungroup() %>%
  distinct(mod, trt_effect_comp, .keep_all = T) %>%
  select(mod, trt_effect_comp, N) %>%
  pivot_wider(names_from = trt_effect_comp, values_from = N)

print("Sample sizes needed for 80% power by model:")
print(sample_sizes_by_model)


# 7. SAVE OUTPUT ===============================================================

# Save the sensitivity analysis plot
png(here("Plots", "R1", "sensitivity_analysis.png"),  width = 15, height = 6, unit = "in", res = 350)
sensitivity_plot
dev.off()

cat("Sensitivity analysis completed!\n")
cat("Plot saved to:", here("Plots", "sensitivity_analysis.png"), "\n")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
