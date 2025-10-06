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
power_threshold <- 0.8

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


all_estimation_error <-  all_summary_stats %>%
  group_by(day_plans, N, N_swabs_per_day, trt_effect_comp, 
           k_slope, k_sigma_logvl, k_sigmasq_u_1, k_sigmasq_u_2) %>%
  mutate(n_per_group = n(),
         sim_k2 = 1:n_per_group[1],
         median_effect = (q2),
         diff_estimates = median_effect - log(trt_effect_comp)) 

all_estimation_error$trt_effect_comp <- as.factor(all_estimation_error$trt_effect_comp)
levels(all_estimation_error$trt_effect_comp) <-  c("Treatment effect = 0%",
                                          "Treatment effect = 20%",
                                          "Treatment effect = 40%",
                                          "Treatment effect = 60%",
                                          "Treatment effect = 80%",
                                          "Treatment effect = 100%")

#  Baseline ==============================================================
plot_error_dat <- all_estimation_error %>% filter(
  k_slope %in% 1,
  k_sigma_logvl %in% 1,
  k_sigmasq_u_2 %in% 1,
  k_sigmasq_u_1 %in% 1,
  N_swabs_per_day %in% 2,
  day_plans %in% "0,1,2,3,4,5",
  N %in% c(10,40,160,640)
) 

med_error <- plot_error_dat %>%
  summarise(med = median(diff_estimates),
            low_IQR = quantile(diff_estimates, 0.25),
            up_IQR = quantile(diff_estimates, 0.75),
            min = min(diff_estimates),
            max = max(diff_estimates))

G1 <- med_error %>%
  ggplot() +
  scale_x_log10(breaks = c(10, 40, 160, 640)) +
  theme_bw(base_size = 15) +
  geom_hline(yintercept = 0, col = 'red',linetype = "22", linewidth = 0.7, alpha = 0.8) +
  scale_color_manual(values = rev(c("#1230AE", "#6C48C5", "#C68FE6","#D8A25E", "#A04747", "black")),  name = "") +
  geom_jitter(aes(x = N, y = med,  col = as.factor(trt_effect_comp)), size = 4,
              position = position_dodge(width = 0.3), alpha = 0.85) +
  geom_errorbar(aes(x = N, ymin = low_IQR, ymax = up_IQR, col = as.factor(trt_effect_comp)),
                width = 0, position = position_dodge(width = 0.3), alpha = 0.75, linewidth = 0.75) +
  xlab("Number of patients per arm") +
  ylab("Estimation error") +
  ggtitle("A) Baseline viral kinetics") +
  theme(axis.title = element_text(face = "bold")) +
  ylim(-0.8, 0.3) 
  
G1
# B) INTERCEPT VARIATION -------------------------------------------------------
plot_error_dat <- all_estimation_error %>% filter(
  k_slope %in% 1,
  k_sigma_logvl %in% 1,
  !k_sigmasq_u_1 %in% 1,
  k_sigmasq_u_2 %in% 1,
  N_swabs_per_day %in% 2,
  day_plans %in% "0,1,2,3,4,5",
  N %in% c(10, 40, 160, 640)
) 

plot_error_dat$k_sigmasq_u_1 <-  as.factor(plot_error_dat$k_sigmasq_u_1)
levels(plot_error_dat$k_sigmasq_u_1) <- c("0.5x", "2.0x")


med_error <- plot_error_dat %>%
  summarise(med = median(diff_estimates),
            low_IQR = quantile(diff_estimates, 0.25),
            up_IQR = quantile(diff_estimates, 0.75),
            min = min(diff_estimates),
            max = max(diff_estimates))

G2<-med_error %>%
  ggplot() +
  scale_x_log10(breaks = c(10, 40, 160, 640)) +
  theme_bw(base_size = 12) +
  geom_hline(yintercept = 0, col = 'red',linetype = "22", linewidth = 0.5, alpha = 0.8) +
  scale_color_manual(values = rev(c("#1230AE", "#6C48C5", "#C68FE6","#D8A25E", "#A04747", "black")),  name = "") +
  geom_jitter(aes(x = N, y = med,  col = as.factor(trt_effect_comp)), size = 1.4,
              position = position_dodge(width = 0.3), alpha = 0.85) +
  geom_errorbar(aes(x = N, ymin = low_IQR, ymax = up_IQR, col = as.factor(trt_effect_comp)),
                width = 0, position = position_dodge(width = 0.3), alpha = 0.75, linewidth = 0.55) +
  xlab("") +
  ylab("") +
  ggtitle(expression("B) Inter-individual variation on the intercept (" * sigma[theta[1]] * ")")) +
  theme(axis.title = element_text(face = "bold")) +
  ylim(-0.8, 0.3) +
  facet_grid(.~k_sigmasq_u_1)
G2
# C) Slope VARIATION -------------------------------------------------------

plot_error_dat <- all_estimation_error %>% filter(
  k_slope %in% 1,
  k_sigma_logvl %in% 1,
  k_sigmasq_u_1 %in% 1,
  !k_sigmasq_u_2 %in% 1,
  N_swabs_per_day %in% 2,
  day_plans %in% "0,1,2,3,4,5",
  N %in% c(10, 40, 160, 640)
) 

plot_error_dat$k_sigmasq_u_2 <-  as.factor(plot_error_dat$k_sigmasq_u_2)
levels(plot_error_dat$k_sigmasq_u_2) <- c("0.5x", "2.0x")


med_error <- plot_error_dat %>%
  summarise(med = median(diff_estimates),
            low_IQR = quantile(diff_estimates, 0.25),
            up_IQR = quantile(diff_estimates, 0.75),
            min = min(diff_estimates),
            max = max(diff_estimates))

G3 <- med_error %>%
  ggplot() +
  scale_x_log10(breaks = c(10, 40, 160, 640)) +
  theme_bw(base_size = 12) +
  geom_hline(yintercept = 0, col = 'red',linetype = "22", linewidth = 0.5, alpha = 0.8) +
  scale_color_manual(values = rev(c("#1230AE", "#6C48C5", "#C68FE6","#D8A25E", "#A04747", "black")),  name = "") +
  geom_jitter(aes(x = N, y = med,  col = as.factor(trt_effect_comp)), size = 1.4,
              position = position_dodge(width = 0.3), alpha = 0.85) +
  geom_errorbar(aes(x = N, ymin = low_IQR, ymax = up_IQR, col = as.factor(trt_effect_comp)),
                width = 0, position = position_dodge(width = 0.3), alpha = 0.75, linewidth = 0.55) +
  xlab("") +
  ylab("") +
  ggtitle(expression("C) Inter-individual variation on the slope (" * sigma[theta[2]] * ")")) +
  theme(axis.title = element_text(face = "bold")) +
  ylim(-0.8, 0.3) +
  facet_grid(.~k_sigmasq_u_2)
G3
# D) Observation VARIATION -------------------------------------------------------

plot_error_dat <- all_estimation_error %>% filter(
  k_slope %in% 1,
  !k_sigma_logvl %in% 1,
  k_sigmasq_u_1 %in% 1,
  k_sigmasq_u_2 %in% 1,
  N_swabs_per_day %in% 2,
  day_plans %in% "0,1,2,3,4,5",
  N %in% c(10, 40, 160, 640)
) 

plot_error_dat$k_sigma_logvl <-  as.factor(plot_error_dat$k_sigma_logvl)
levels(plot_error_dat$k_sigma_logvl) <- c("0.5x", "2.0x")


med_error <- plot_error_dat %>%
  summarise(med = median(diff_estimates),
            low_IQR = quantile(diff_estimates, 0.25),
            up_IQR = quantile(diff_estimates, 0.75),
            min = min(diff_estimates),
            max = max(diff_estimates))

G4 <- med_error %>%
  ggplot() +
  scale_x_log10(breaks = c(10, 40, 160, 640)) +
  theme_bw(base_size = 12) +
  geom_hline(yintercept = 0, col = 'red',linetype = "22", linewidth = 0.5, alpha = 0.8) +
  scale_color_manual(values = rev(c("#1230AE", "#6C48C5", "#C68FE6","#D8A25E", "#A04747", "black")),  name = "") +
  geom_jitter(aes(x = N, y = med,  col = as.factor(trt_effect_comp)), size = 1.4,
              position = position_dodge(width = 0.3), alpha = 0.85) +
  geom_errorbar(aes(x = N, ymin = low_IQR, ymax = up_IQR, col = as.factor(trt_effect_comp)),
                width = 0, position = position_dodge(width = 0.3), alpha = 0.75, linewidth = 0.55) +
  xlab("") +
  ylab("") +
  ggtitle(expression("D) Observation variation (" * sigma[VL] * ")")) +
  theme(axis.title = element_text(face = "bold")) +
  ylim(-0.8, 0.3) +
  facet_grid(.~k_sigma_logvl)
G4

# E) SAMPLING SCHEDULE ---------------------------------------------------------
plot_error_dat <- all_estimation_error %>% filter(
  k_slope %in% 1,
  k_sigma_logvl %in% 1,
  k_sigmasq_u_1 %in% 1,
  k_sigmasq_u_2 %in% 1,
  N_swabs_per_day %in% 2,
  #day_plans %in% "0,1,2,3,4,5",
  N %in% c(10, 40, 160, 640)
) 

plot_error_dat$day_plans <-  as.factor(plot_error_dat$day_plans)
levels(plot_error_dat$day_plans) <- c("Everyday\n(n=6)", "Every other day\n(n=3)", "First and last day\n(n=2)")

med_error <- plot_error_dat %>%
  summarise(med = median(diff_estimates),
            low_IQR = quantile(diff_estimates, 0.25),
            up_IQR = quantile(diff_estimates, 0.75),
            min = min(diff_estimates),
            max = max(diff_estimates))

G5 <- med_error %>%
  ggplot() +
  scale_x_log10(breaks = c(10, 40, 160, 640)) +
  theme_bw(base_size = 12) +
  geom_hline(yintercept = 0, col = 'red',linetype = "22", linewidth = 0.5, alpha = 0.8) +
  scale_color_manual(values = rev(c("#1230AE", "#6C48C5", "#C68FE6","#D8A25E", "#A04747", "black")),  name = "") +
  geom_jitter(aes(x = N, y = med,  col = as.factor(trt_effect_comp)), size = 1.4,
              position = position_dodge(width = 0.3), alpha = 0.85) +
  geom_errorbar(aes(x = N, ymin = low_IQR, ymax = up_IQR, col = as.factor(trt_effect_comp)),
                width = 0, position = position_dodge(width = 0.3), alpha = 0.75, linewidth = 0.55) +
  xlab("") +
  ylab("") +
  ggtitle("A) Varied sampling schedule") +
  theme(axis.title = element_text(face = "bold")) +
  ylim(-0.8, 0.3) +
  facet_grid(.~day_plans)
G5


# F) NUMBER OF SWABS PER DAY ---------------------------------------------------
plot_error_dat <- all_estimation_error %>% filter(
  k_slope %in% 1,
  k_sigma_logvl %in% 1,
  k_sigmasq_u_1 %in% 1,
  k_sigmasq_u_2 %in% 1,
  #N_swabs_per_day %in% 2,
  day_plans %in% "0,1,2,3,4,5",
  N %in% c(10, 40, 160, 640)
) 

plot_error_dat$N_swabs_per_day <-  as.factor(plot_error_dat$N_swabs_per_day)
levels(plot_error_dat$N_swabs_per_day) <- c("1 swab/day\n", "2 swabs/day\n",  "4 swabs/day\n")

med_error <- plot_error_dat %>%
  summarise(med = median(diff_estimates),
            low_IQR = quantile(diff_estimates, 0.25),
            up_IQR = quantile(diff_estimates, 0.75),
            min = min(diff_estimates),
            max = max(diff_estimates))

G6 <- med_error %>%
  ggplot() +
  scale_x_log10(breaks = c(10, 40, 160, 640)) +
  theme_bw(base_size = 12) +
  geom_hline(yintercept = 0, col = 'red',linetype = "22", linewidth = 0.5, alpha = 0.8) +
  scale_color_manual(values = rev(c("#1230AE", "#6C48C5", "#C68FE6","#D8A25E", "#A04747", "black")),  name = "") +
  geom_jitter(aes(x = N, y = med,  col = as.factor(trt_effect_comp)), size = 1.4,
              position = position_dodge(width = 0.3), alpha = 0.85) +
  geom_errorbar(aes(x = N, ymin = low_IQR, ymax = up_IQR, col = as.factor(trt_effect_comp)),
                width = 0, position = position_dodge(width = 0.3), alpha = 0.75, linewidth = 0.55) +
  xlab("") +
  ylab("") +
  ggtitle("B) Varied number of swabs per day") +
  theme(axis.title = element_text(face = "bold")) +
  ylim(-0.8, 0.3) +
  facet_grid(.~N_swabs_per_day)
G6


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
    bottom = textGrob("Number of patients per arm", vjust = 0.5, gp = gpar(cex = 1.2, fontface="bold")),
    left = textGrob("Estimation error", rot = 90, gp = gpar(cex = 1.2, fontface="bold"))
  )

# Combine baseline with variation plots
combined_plot_all <- ggarrange(
  G1 + theme(legend.position = "none", axis.title = element_text(face = "bold")), 
  combined_plot
)

# Save main comparison plot
png(here("Plots", "R1", "variations_vs_estimation_error.png"), width = 14, height = 8, units = 'in', res = 350)
print(combined_plot_all)
dev.off()


# Combine sampling plots (E, F)
combined_plot2 <- ggarrange(G5, G6, 
                            common.legend = T, ncol = 1,
                            legend = "right") %>%
  annotate_figure(
    bottom = textGrob("Number of patients per arm", vjust = 0.5, gp = gpar(cex = 1.2, fontface="bold")),
    left = textGrob("Estimation error", rot = 90, gp = gpar(cex = 1.2, fontface="bold"))
  )

# Save sampling comparison plot
png(here("Plots", "R1", "sampling_vs_estimation_error.png"), width = 10, height = 8, unit = "in", res = 350)
combined_plot2
dev.off()



