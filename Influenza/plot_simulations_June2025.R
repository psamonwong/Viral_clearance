# ==============================================================================
# POWER ANALYSIS VISUALIZATION FOR VIRAL CLEARANCE STUDY
# ==============================================================================
# This script generates power curves showing probability of detecting treatment 
# effects across different sample sizes and simulation scenarios
# ==============================================================================

# 1. SETUP & DEPENDENCIES =====================================================
# Load required libraries
library(tidyverse)    # Data manipulation and plotting
library(here)         # File path management
library(mgcv)         # GAM modeling
library(ggpubr)       # Plot arrangement
library(grid)         # Grid graphics for annotations

# 2. CONFIGURATION ============================================================
# Set quantiles used for summarizing simulation results
QUANTILES <- c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
rerun <- FALSE
power_threshold <- 0.9

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

# Prepare data for GAM modeling (add power detection flag)
prepare_gam_data <- function(df_summary_all) {
  df_summary_all %>%
    group_by(day_plans, N, N_swabs_per_day, trt_effect_comp, 
             k_slope, k_sigma_logvl, k_sigmasq_u_1, k_sigmasq_u_2
             ) %>%
    mutate(rejectH0 = low95 > 0) 
}

# 4. DATA LOADING ==============================================================
# Load negative control simulation data
all_summary_stats <- load_simulation(
  rdata = 'sim_settings_Jun25.RData',
  simdir = 'sims_out_June2025',
  neg_control = TRUE,
  rerun = F,
  output_dir = "sim_posterior_Jun25_negative_control.csv"
) 

# Prepare data for power analysis
gam_data <- prepare_gam_data(all_summary_stats)

# 5. POWER ANALYSIS PLOTS ======================================================

# A) BASELINE SCENARIO --------------------------------------------------------
# Shows power curves under standard conditions (all parameters at baseline)
dat <- gam_data
dat_for_fit <- dat %>% filter(
  k_slope %in% 1,
  k_sigma_logvl %in% 1,
  k_sigmasq_u_1 %in% 1,
  k_sigmasq_u_2 %in% 1,
  N_swabs_per_day %in% 2,
  day_plans %in% "0,1,2,3,4,5"
)

# Format treatment effect labels for plotting
dat_for_fit$trt_effect_comp <- as.factor(dat_for_fit$trt_effect_comp)
levels(dat_for_fit$trt_effect_comp) <- c(
  "Treatment effect = 0%", "Treatment effect = 20%",
  "Treatment effect = 40%", "Treatment effect = 60%", 
  "Treatment effect = 80%", "Treatment effect = 100%"
)

# Calculate empirical power at each sample size
summary_p <- dat_for_fit %>%
  group_by(N, trt_effect_comp) %>%
  summarise(n = n(), p_rejectH0 = sum(rejectH0) / n)

# Fit GAM model for smooth power curves
fit <- mgcv::gam(rejectH0 ~ log10(N)*trt_effect_comp, data = dat_for_fit, family = "binomial")
new_dat <- expand.grid(
  N = seq(1,2000,1),
  trt_effect_comp = levels(as.factor(dat_for_fit$trt_effect_comp))
)
new_dat$preds = predict(fit, newdata = new_dat, type = 'response') %>% as.numeric()

# Create baseline power plot
G1 <- ggplot(mapping = aes(col = trt_effect_comp)) +
  geom_line(data = new_dat, mapping = aes(x = N, y = preds),
            linewidth = 1, alpha = 0.75, linetype = '31') +
  geom_point(data = summary_p, mapping = aes(x = N, y = p_rejectH0),
             size = 3, alpha = 0.75) +
  scale_x_log10(breaks = c(10, 20, 40, 80, 160, 320, 640), limits = c(10, 640)) +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  theme_bw(base_size = 15) +
  scale_color_manual(values = rev(c("#1230AE", "#6C48C5", "#C68FE6","#D8A25E", "#A04747", "black")), name = "") +
  geom_hline(yintercept = c(0.8, 0.9), col = "gray20", linetype = "22", linewidth = 0.75, alpha = 0.8) +
  geom_hline(yintercept = c(0,1), col = "grey20", linetype = "22", linewidth = 0.75, alpha = 0.2) +
  xlab("Number of patients per arm") +
  ylab("Pr[Treatment effect ≤ 0] < 0.025") +
  ggtitle("A) Baseline viral kinetics")
G1
# Calculate sample sizes needed for 80% power
new_dat_diff <- new_dat %>%
  mutate(diff_pred = abs(preds-power_threshold)) %>%
  group_by(trt_effect_comp) %>%
  filter(diff_pred == min(diff_pred) & diff_pred < 0.01) %>%
  ungroup() %>%
  print()

# B) INTERCEPT VARIATION -------------------------------------------------------
# Compare power under different levels of inter-individual variation in intercept
dat <- gam_data
dat_for_fit <- dat %>% filter(
  k_slope %in% 1,
  k_sigma_logvl %in% 1,
  k_sigmasq_u_2 %in% 1,
  N_swabs_per_day %in% 2,
  day_plans %in% "0,1,2,3,4,5"
)

dat_for_fit$trt_effect_comp <- as.factor(dat_for_fit$trt_effect_comp)
levels(dat_for_fit$trt_effect_comp) <-  c("Treatment effect = 0%",
                                          "Treatment effect = 20%",
                                          "Treatment effect = 40%",
                                          "Treatment effect = 60%",
                                          "Treatment effect = 80%",
                                          "Treatment effect = 100%")
dat_for_fit$k_sigmasq_u_1 <-  as.factor(dat_for_fit$k_sigmasq_u_1)
levels(dat_for_fit$k_sigmasq_u_1) <- c("0.5x", "1.0x", "2.0x")

summary_p <- dat_for_fit %>%
  group_by(N, trt_effect_comp, k_sigmasq_u_1) %>%
  summarise(
    n = n(),
    p_rejectH0 = sum(rejectH0) / n)

fit <- mgcv::gam(rejectH0 ~ log10(N)*trt_effect_comp*k_sigmasq_u_1,  data = dat_for_fit, family = "binomial")
new_dat <- expand.grid(N = seq(1,2000,1),
                       k_sigmasq_u_1 = levels(as.factor(dat_for_fit$k_sigmasq_u_1)),
                       trt_effect_comp =  levels(as.factor(dat_for_fit$trt_effect_comp)))
new_dat$preds = predict(fit, newdata = new_dat, type = 'response') %>% as.numeric()

new_dat <- new_dat %>% filter(k_sigmasq_u_1 != "1.0x")
summary_p <- summary_p %>% filter(k_sigmasq_u_1 != "1.0x")


G2 <- ggplot(mapping = aes(col = trt_effect_comp)) +
  geom_line(data = new_dat, mapping = aes(x = N, y = preds),
            linewidth = 0.75, alpha = 0.75, linetype = '31') +
  geom_point(data = summary_p, mapping = aes(x = N, y = p_rejectH0),
             size = 1.4, alpha = 0.75) +
  scale_x_log10(breaks = c(10, 20, 40, 80, 160, 320, 640), limits = c(10, 640)) +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  facet_wrap(.~k_sigmasq_u_1, nrow = 1) +
  theme_bw(base_size = 11) +
  scale_color_manual(values = rev(c("#1230AE", "#6C48C5", "#C68FE6","#D8A25E", "#A04747", "black")),  name = "") +
  geom_hline(yintercept = c(0.8, 0.9), col = "gray20", linetype = "22", linewidth = 0.75, alpha = 0.8) +
  geom_hline(yintercept = c(0,1), col = "grey20", linetype = "22", linewidth = 0.4, alpha = 0.2) +
  xlab("") +
  ylab("") +
  ggtitle(expression("B) Inter-individual variation on the intercept (" * sigma[theta[1]] * ")")) 
G2   

new_dat_diff <- new_dat %>%
  mutate(diff_pred = abs(preds-power_threshold)) %>%
  group_by(trt_effect_comp,k_sigmasq_u_1) %>%
  filter(diff_pred == min(diff_pred) &
           diff_pred < 0.01) %>%
  ungroup() %>%
  distinct(k_sigmasq_u_1, trt_effect_comp, .keep_all = T) %>%
  select(k_sigmasq_u_1, trt_effect_comp, N) %>%
  pivot_wider(names_from = trt_effect_comp, values_from = N) %>%
  print() 


# C) SLOPE VARIATION -----------------------------------------------------------
# Compare power under different levels of inter-individual variation in slope
dat <- gam_data
dat_for_fit <- dat %>% filter(
  k_slope %in% 1,
  k_sigma_logvl %in% 1,
  k_sigmasq_u_1 %in% 1,
  N_swabs_per_day %in% 2,
  day_plans %in% "0,1,2,3,4,5"
)

dat_for_fit$trt_effect_comp <- as.factor(dat_for_fit$trt_effect_comp)
levels(dat_for_fit$trt_effect_comp) <-  c("Treatment effect = 0%",
                                          "Treatment effect = 20%",
                                          "Treatment effect = 40%",
                                          "Treatment effect = 60%",
                                          "Treatment effect = 80%",
                                          "Treatment effect = 100%")
dat_for_fit$k_sigmasq_u_2 <-  as.factor(dat_for_fit$k_sigmasq_u_2)
levels(dat_for_fit$k_sigmasq_u_2) <- c("0.5x", "1.0x", "2.0x")

summary_p <- dat_for_fit %>%
  group_by(N, trt_effect_comp, k_sigmasq_u_2) %>%
  summarise(
    n = n(),
    p_rejectH0 = sum(rejectH0) / n)

fit <- mgcv::gam(rejectH0 ~ log10(N)*trt_effect_comp*k_sigmasq_u_2,  data = dat_for_fit, family = "binomial")
new_dat <- expand.grid(N = seq(1,2000,1),
                       k_sigmasq_u_2 = levels(as.factor(dat_for_fit$k_sigmasq_u_2)),
                       trt_effect_comp =  levels(as.factor(dat_for_fit$trt_effect_comp)))
new_dat$preds = predict(fit, newdata = new_dat, type = 'response') %>% as.numeric()

new_dat <- new_dat %>% filter(k_sigmasq_u_2 != "1.0x")
summary_p <- summary_p %>% filter(k_sigmasq_u_2 != "1.0x")

G3 <- ggplot(mapping = aes(col = trt_effect_comp)) +
  geom_line(data = new_dat, mapping = aes(x = N, y = preds),
            linewidth = 0.75, alpha = 0.75, linetype = '31') +
  geom_point(data = summary_p, mapping = aes(x = N, y = p_rejectH0),
             size = 1.4, alpha = 0.75) +
  scale_x_log10(breaks = c(10, 20, 40, 80, 160, 320, 640), limits = c(10, 640)) +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  facet_wrap(.~k_sigmasq_u_2, nrow = 1) +
  theme_bw(base_size = 11) +
  scale_color_manual(values = rev(c("#1230AE", "#6C48C5", "#C68FE6","#D8A25E", "#A04747", "black")),  name = "") +
  geom_hline(yintercept = c(0.8, 0.9), col = "gray20", linetype = "22", linewidth = 0.75, alpha = 0.8) +
  geom_hline(yintercept = c(0,1), col = "grey20", linetype = "22", linewidth = 0.4, alpha = 0.2) +
  xlab("") +
  ylab("") +
  ggtitle(expression("C) Inter-individual variation on the slope (" * sigma[theta[2]] * ")")) 
G3

new_dat_diff <- new_dat %>%
  mutate(diff_pred = abs(preds-power_threshold)) %>%
  group_by(trt_effect_comp,k_sigmasq_u_2) %>%
  filter(diff_pred == min(diff_pred) &
           diff_pred < 0.01) %>%
  ungroup() %>%
  distinct(k_sigmasq_u_2, trt_effect_comp, .keep_all = T) %>%
  select(k_sigmasq_u_2, trt_effect_comp, N) %>%
  pivot_wider(names_from = trt_effect_comp, values_from = N) %>%
  print()

# D) OBSERVATIONAL ERROR -------------------------------------------------------
# Compare power under different levels of measurement error
dat <- gam_data
dat_for_fit <- dat %>% filter(
  k_slope %in% 1,
  k_sigmasq_u_1 %in% 1,
  k_sigmasq_u_2 %in% 1,
  N_swabs_per_day %in% 2,
  day_plans %in% "0,1,2,3,4,5"
)

dat_for_fit$trt_effect_comp <- as.factor(dat_for_fit$trt_effect_comp)
levels(dat_for_fit$trt_effect_comp) <-  c("Treatment effect = 0%",
                                          "Treatment effect = 20%",
                                          "Treatment effect = 40%",
                                          "Treatment effect = 60%",
                                          "Treatment effect = 80%",
                                          "Treatment effect = 100%")
dat_for_fit$k_sigma_logvl <-  as.factor(dat_for_fit$k_sigma_logvl)
levels(dat_for_fit$k_sigma_logvl) <- c("0.5x", "1.0x", "2.0x")

summary_p <- dat_for_fit %>%
  group_by(N, trt_effect_comp, k_sigma_logvl) %>%
  summarise(
    n = n(),
    p_rejectH0 = sum(rejectH0) / n)

fit <- mgcv::gam(rejectH0 ~ log10(N)*trt_effect_comp*k_sigma_logvl,  data = dat_for_fit, family = "binomial")
new_dat <- expand.grid(N = seq(1,2000,1),
                       k_sigma_logvl = levels(as.factor(dat_for_fit$k_sigma_logvl)),
                       trt_effect_comp =  levels(as.factor(dat_for_fit$trt_effect_comp)))
new_dat$preds = predict(fit, newdata = new_dat, type = 'response') %>% as.numeric()

new_dat <- new_dat %>% filter(k_sigma_logvl != "1.0x")
summary_p <- summary_p %>% filter(k_sigma_logvl != "1.0x")

G4 <- ggplot(mapping = aes(col = trt_effect_comp)) +
  geom_line(data = new_dat, mapping = aes(x = N, y = preds),
            linewidth = 0.75, alpha = 0.75, linetype = '31') +
  geom_point(data = summary_p, mapping = aes(x = N, y = p_rejectH0),
             size = 1.4, alpha = 0.75) +
  scale_x_log10(breaks = c(10, 20, 40, 80, 160, 320, 640), limits = c(10, 640)) +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  facet_wrap(.~k_sigma_logvl, nrow = 1) +
  theme_bw(base_size = 11) +
  scale_color_manual(values = rev(c("#1230AE", "#6C48C5", "#C68FE6","#D8A25E", "#A04747", "black")),  name = "") +
  geom_hline(yintercept = c(0.8, 0.9), col = "gray20", linetype = "22", linewidth = 0.75, alpha = 0.8) +
  geom_hline(yintercept = c(0,1), col = "grey20", linetype = "22", linewidth = 0.4, alpha = 0.2) +
  xlab("") +
  ylab("") +
  ggtitle(expression("D) Observation variation (" * sigma[VL] * ")")) 
G4

new_dat_diff <- new_dat %>%
  mutate(diff_pred = abs(preds-power_threshold)) %>%
  group_by(trt_effect_comp,k_sigma_logvl) %>%
  filter(diff_pred == min(diff_pred) &
           diff_pred < 0.01) %>%
  ungroup() %>%
  distinct(k_sigma_logvl, trt_effect_comp, .keep_all = T) %>%
  select(k_sigma_logvl, trt_effect_comp, N) %>%
  pivot_wider(names_from = trt_effect_comp, values_from = N) %>%
  print()
  
# E) SAMPLING SCHEDULE ---------------------------------------------------------
# Compare power across different sampling schedules
dat <- gam_data
dat_for_fit <- dat %>% filter(
  k_slope %in% 1,
  k_sigma_logvl %in% 1,
  k_sigmasq_u_1 %in% 1,
  k_sigmasq_u_2 %in% 1,
  N_swabs_per_day %in% 2
)

dat_for_fit$trt_effect_comp <- as.factor(dat_for_fit$trt_effect_comp)
levels(dat_for_fit$trt_effect_comp) <-  c("Treatment effect = 0%",
                                          "Treatment effect = 20%",
                                          "Treatment effect = 40%",
                                          "Treatment effect = 60%",
                                          "Treatment effect = 80%",
                                          "Treatment effect = 100%")
dat_for_fit$day_plans <-  as.factor(dat_for_fit$day_plans)
levels(dat_for_fit$day_plans) <- c("Everyday\n(n=6)", "Every other day\n(n=3)", "First and last day\n(n=2)")

summary_p <- dat_for_fit %>%
  group_by(N, trt_effect_comp, day_plans) %>%
  summarise(
    n = n(),
    p_rejectH0 = sum(rejectH0) / n)

fit <- mgcv::gam(rejectH0 ~ log10(N)*trt_effect_comp*day_plans,  data = dat_for_fit, family = "binomial")
new_dat <- expand.grid(N = seq(1,2000,1),
                       day_plans = levels(as.factor(dat_for_fit$day_plans)),
                       trt_effect_comp =  levels(as.factor(dat_for_fit$trt_effect_comp)))
new_dat$preds = predict(fit, newdata = new_dat, type = 'response') %>% as.numeric()

G5 <- ggplot(mapping = aes(col = trt_effect_comp)) +
  geom_line(data = new_dat, mapping = aes(x = N, y = preds),
            linewidth = 0.8, alpha = 0.75, linetype = '31') +
  geom_point(data = summary_p, mapping = aes(x = N, y = p_rejectH0),
             size = 1.8, alpha = 0.75) +
  scale_x_log10(breaks = c(10, 20, 40, 80, 160, 320, 640), limits = c(10, 640)) +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  facet_wrap(.~day_plans, nrow = 1) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(c("#1230AE", "#6C48C5", "#C68FE6","#D8A25E", "#A04747", "black")),  name = "") +
  geom_hline(yintercept = c(0.8, 0.9), col = "gray20", linetype = "22", linewidth = 0.75, alpha = 0.8) +
  geom_hline(yintercept = c(0,1), col = "grey20", linetype = "22", linewidth = 0.6, alpha = 0.2) +
  xlab("") +
  ylab("") +
  ggtitle("A) Varied sampling schedule") 
G5 
new_dat_diff <- new_dat %>%
  mutate(diff_pred = abs(preds-power_threshold)) %>%
  group_by(trt_effect_comp,day_plans) %>%
  filter(diff_pred == min(diff_pred) &
           diff_pred < 0.01) %>%
  ungroup() %>%
  distinct(day_plans, trt_effect_comp, .keep_all = T) %>%
  select(day_plans, trt_effect_comp, N) %>%
  pivot_wider(names_from = trt_effect_comp, values_from = N) %>%
  print()


# F) NUMBER OF SWABS PER DAY ---------------------------------------------------
# Compare power with different numbers of daily measurements
dat <- gam_data
dat_for_fit <- dat %>% filter(
  k_slope %in% 1,
  k_sigma_logvl %in% 1,
  k_sigmasq_u_1 %in% 1,
  k_sigmasq_u_2 %in% 1,
  day_plans %in% "0,1,2,3,4,5"
)

dat_for_fit$trt_effect_comp <- as.factor(dat_for_fit$trt_effect_comp)
levels(dat_for_fit$trt_effect_comp) <-  c("Treatment effect = 0%",
                                          "Treatment effect = 20%",
                                          "Treatment effect = 40%",
                                          "Treatment effect = 60%",
                                          "Treatment effect = 80%",
                                          "Treatment effect = 100%")
dat_for_fit$N_swabs_per_day <-  as.factor(dat_for_fit$N_swabs_per_day)
levels(dat_for_fit$N_swabs_per_day) <- c("1 swab/day\n", "2 swabs/day\n",  "4 swabs/day\n")

summary_p <- dat_for_fit %>%
  group_by(N, trt_effect_comp, N_swabs_per_day) %>%
  summarise(
    n = n(),
    p_rejectH0 = sum(rejectH0) / n)

fit <- mgcv::gam(rejectH0 ~ log10(N)*trt_effect_comp*N_swabs_per_day,  data = dat_for_fit, family = "binomial")
new_dat <- expand.grid(N = seq(1,2000,1),
                       N_swabs_per_day = levels(as.factor(dat_for_fit$N_swabs_per_day)),
                       trt_effect_comp =  levels(as.factor(dat_for_fit$trt_effect_comp)))
new_dat$preds = predict(fit, newdata = new_dat, type = 'response') %>% as.numeric()

G6 <- ggplot(mapping = aes(col = trt_effect_comp)) +
  geom_line(data = new_dat, mapping = aes(x = N, y = preds),
            linewidth = 0.8, alpha = 0.75, linetype = '31') +
  geom_point(data = summary_p, mapping = aes(x = N, y = p_rejectH0),
             size = 1.8, alpha = 0.75) +
  scale_x_log10(breaks = c(10, 20, 40, 80, 160, 320, 640), limits = c(10, 640)) +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  facet_wrap(.~N_swabs_per_day, nrow = 1) +
  theme_bw(base_size = 12) +
  scale_color_manual(values = rev(c("#1230AE", "#6C48C5", "#C68FE6","#D8A25E", "#A04747", "black")),  name = "") +
  geom_hline(yintercept = c(0.8, 0.9), col = "gray20", linetype = "22", linewidth = 0.75, alpha = 0.8) +
  geom_hline(yintercept = c(0,1), col = "grey20", linetype = "22", linewidth = 0.6, alpha = 0.2) +
  xlab("") +
  ylab("") +
  ggtitle("B) Varied number of swabs per day") 
G6

new_dat_diff <- new_dat %>%
  mutate(diff_pred = abs(preds-power_threshold)) %>%
  group_by(trt_effect_comp,N_swabs_per_day) %>%
  filter(diff_pred == min(diff_pred) &
           diff_pred < 0.01) %>%
  ungroup() %>%
  distinct(N_swabs_per_day, trt_effect_comp, .keep_all = T) %>%
  select(N_swabs_per_day, trt_effect_comp, N) %>%
  pivot_wider(names_from = trt_effect_comp, values_from = N) %>%
  print()



# 6. PLOT ASSEMBLY & OUTPUT ====================================================

# Prepare plots for combination (adjust margins)
G2 <- G2 + theme(plot.margin = margin(5, 5, 5, 5, "pt"))
G3 <- G3 + theme(plot.margin = margin(5, 5, 5, 5, "pt"))
G4 <- G4 + theme(plot.margin = margin(5, 5, 5, 5, "pt"))

# Combine variation plots (B, C, D)
combined_plot <- ggarrange(G2, G3, G4, 
                           common.legend = T, ncol = 1,
                           legend = "right", align = "hv") %>%
  annotate_figure(
    bottom = textGrob("Number of patients per arm", vjust = 0.5, gp = gpar(cex = 1.2)),
    left = textGrob("Pr[Treatment effect ≤ 0] < 0.025", rot = 90, gp = gpar(cex = 1.2))
  )

# Combine baseline with variation plots
combined_plot_all <- ggarrange(
  G1 + theme(legend.position = "none"), 
  combined_plot
)

# Save main comparison plot
png(here("Plots", "R1", "variations_vs_power.png"), width = 14, height = 8, units = 'in', res = 350)
print(combined_plot_all)
dev.off()

# Prepare sampling plots
G5 <- G5 + theme(plot.margin = margin(5, 5, 5, 5, "pt"))
G6 <- G6 + theme(plot.margin = margin(5, 5, 5, 5, "pt"))

# Combine sampling plots (E, F)
combined_plot2 <- ggarrange(G5, G6, 
                           common.legend = T, ncol = 1,
                           legend = "right") %>%
  annotate_figure(
    bottom = textGrob("Number of patients per arm", vjust = 0.5, gp = gpar(cex = 1.2)),
    left = textGrob("Pr[Treatment effect ≤ 0] < 0.025", rot = 90, gp = gpar(cex = 1.2))
  )

# Save sampling comparison plot
png(here("Plots", "R1", "sampling_vs_power.png"), width = 10, height = 8, unit = "in", res = 350)
combined_plot2
dev.off()

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
