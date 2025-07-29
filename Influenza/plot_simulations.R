# 1. Load required libraries
# Load core tidyverse for data manipulation and plotting, and here for file path management
library(tidyverse)
library(here)

# 2. Define constants
# Set quantiles used for summarizing simulation results
QUANTILES <- c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
rerun <- FALSE

# 3. Define helper functions
# Relabel factors for consistent and informative plotting
relabel_factors <- function(df) {
  df$day_plans <- as.factor(df$day_plans)
  levels(df$day_plans) <- c("PLATCOV:\n (D0 to D5)", "Mateja et al:\n (D0 and D3)")
  df$N_swabs_per_day <- as.factor(df$N_swabs_per_day)
  levels(df$N_swabs_per_day) <- c("1 swab/day", "2 swabs/day")
  df$trt_effect_comp <- as.factor(df$trt_effect_comp)
  levels(df$trt_effect_comp) <- c(
    "Treatment effect = 0%",
    "Treatment effect = 20%",
    "Treatment effect = 40%",
    "Treatment effect = 60%",
    "Treatment effect = 80%",
    "Treatment effect = 100%"
  )
  df
}

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

# Assign highlight type for plotting emphasis
assign_highlight_type <- function(df) {
  df %>% mutate(type = "Normal",
                type = if_else((as.numeric(day_plans) == 1 & as.numeric(N_swabs_per_day) == 2)|
                                 (as.numeric(day_plans) == 2 & as.numeric(N_swabs_per_day) == 1),
                               "Highlighted", type))
}

# Fit a GAM model and predict power curves
fit_gam_and_predict <- function(dat) {
  fit <- mgcv::gam(rejectH0 ~ log(N)*day_plans*trt_effect_comp*N_swabs_per_day, data = dat, family = "binomial")
  new_dat <- expand.grid(N = seq(10,2000,1),
                        day_plans = levels(dat$day_plans),
                        trt_effect_comp =  levels(dat$trt_effect_comp),
                        N_swabs_per_day =  levels(dat$N_swabs_per_day))
  new_dat <- assign_highlight_type(new_dat)
  new_dat$preds = predict(fit, newdata = new_dat, type = 'response') %>% as.numeric()
  list(fit = fit, new_dat = new_dat)
}

# Calculate power table for each scenario
calculate_power_table <- function(df_summary_all) {
  df_summary_all %>%
    group_by(day_plans, N, N_swabs_per_day, trt_effect_comp) %>%
    mutate(rejectH0 = low95 > 0) %>%
    summarise(power = mean(rejectH0), n = n(), .groups = 'drop') %>%
    relabel_factors() %>%
    assign_highlight_type()
}

# Prepare data for GAM fitting
prepare_gam_data <- function(df_summary_all) {
  df_summary_all %>%
    group_by(day_plans, N, N_swabs_per_day, trt_effect_comp) %>%
    mutate(rejectH0 = low95 > 0) %>%
    relabel_factors()
}

# Summarize estimation error for each scenario
summarize_error_table <- function(df_summary_all) {
  df_summary_all %>%
    group_by(day_plans, N, N_swabs_per_day, trt_effect_comp) %>%
    summarise(med_error = quantile(error, 0.5),
              q1_error = quantile(error, 0.25),
              q3_error = quantile(error, 0.75), .groups = 'drop') %>%
    relabel_factors()
}

# 4. Load data
# Load simulation settings and simulation results
load_simulation <- function(rdata, simdir, rerun, output_dir){
  load(here('Rout', rdata))
  sim_dir <- here(simdir)
  sim_settings <- sim_settings %>% mutate(i = row_number())
  sim_list <- list.files(sim_dir, full.names = TRUE)
  sim_ids <- str_extract(sim_list, "[0-9]+") %>% as.numeric()
  
  if(nrow(sim_settings) == 4600){
    sim_settings <- sim_settings %>%
      filter(i %in%sim_ids)
  }

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


all_summary_stats1 <- load_simulation(rdata = 'sim_settings_sampling_schedule.RData',
                                      simdir = 'sims_out_sampling_schedule',
                                      rerun = FALSE,
                                      output_dir = "sim_posterior.csv") 

all_summary_stats2 <- load_simulation(rdata = 'sim_settings_sampling_schedule_supplement.RData',
                                      simdir = 'sims_out_sampling_schedule_supplement',
                                      rerun = T,
                                      output_dir = "sim_posterior_supplement.csv") 
all_summary_stats <- rbind(all_summary_stats1, all_summary_stats2)

# 6. Calculate power and prepare data for modeling
# Calculate power for each scenario and prepare data for GAM fitting
power_summary <- calculate_power_table(all_summary_stats)
gam_data <- prepare_gam_data(all_summary_stats)
gam_fit_results <- fit_gam_and_predict(gam_data)
predicted_power <- gam_fit_results$new_dat
error_summary <- summarize_error_table(all_summary_stats)

# 7. Plotting and output
# Plot power curve with observed and fitted values
G3 <- power_summary %>%
  filter(!trt_effect_comp %in% c("Treatment effect = 0%", "Treatment effect = 100%")) %>%
  ggplot(aes(x = N, y = power,
             linetype = N_swabs_per_day,
             shape = N_swabs_per_day,
             colour = day_plans,
             alpha = type,
             linewidth = type,
             size = type)) +
  geom_point() +
  scale_alpha_manual(values = c(0.8, 0.3), guide = "none") +
  scale_linewidth_manual(values = c(1, 0.7), guide = "none") +
  scale_size_manual(values = c(3, 2), guide = "none") +
  facet_wrap(.~trt_effect_comp, ncol = 2) +
  theme_bw(base_size = 15) +
  geom_hline(yintercept = c(0.8, 0.9), linetype = '21', col = 'black' , linewidth = 0.7, alpha = 0.5) +
  scale_color_manual(values = (c("#3D90D7", "#B03052")), name = "") +
  scale_linetype_manual(values = c("31", "solid"), name = "") +
  scale_shape_manual(values = c(16,17), name = "") +
  ylab("Probability[Treatment effect > 0] > 0.975") +
  xlab("Number of patients per arm") +
  scale_x_log10(breaks = c(10,  20 , 40, 100, 200 , 400 ,800), limits = c(10,800)) +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  geom_line(data = predicted_power %>%
              filter(!trt_effect_comp %in% c("Treatment effect = 0%", "Treatment effect = 100%")), 
            mapping = aes(x = N, y = preds)) +
  theme(axis.title = element_text(face = 'bold')) +
  guides(shape = guide_legend(override.aes = list(size = 2.5, linewidth = 0.7, alpha = 0.8)),
         color = guide_legend(override.aes = list(size = 2.5, linewidth = 0.7, alpha = 0.8)))
G3
ggsave(here("Plots", "simulation_power_fit.png"), G3, width = 9, height = 8, units = 'in', dpi = 300)

# Find minimum N for each scenario to achieve power threshold
predicting_power <- function(predicted_power, power_threshold){
  levels(predicted_power$day_plans) <- c("PLATCOV", "Mateja")
  levels(predicted_power$N_swabs_per_day) <- c("1swab", "2swabs")
  levels(predicted_power$trt_effect_comp) <- c("0%", "20%","40%","60%","80%","100%")
  res <- predicted_power %>%
    filter(!trt_effect_comp %in% c("0%", "100%")) %>%
    mutate(diff_pred = abs(preds - power_threshold)) %>%
    group_by(day_plans, N_swabs_per_day, trt_effect_comp) %>%
    filter(diff_pred == min(diff_pred) & diff_pred < 0.01) %>%
    ungroup() %>%
    distinct(day_plans, N_swabs_per_day, trt_effect_comp, .keep_all = T) %>%
    select(day_plans, N_swabs_per_day, trt_effect_comp, N) %>%
    pivot_wider(names_from = c(day_plans, N_swabs_per_day), values_from = N)
  res <- res[,c(1,2,4,3,5)]
  colnames(res) <- c("Treatment effect size", 
                     "PLATCOV: 1 swab/day", "PLATCOV: 2 swabs/day",
                     "Mateja: 1 swab/day", "Mateja: 2 swabs/day"
                     )
  res
  }

predicting_power(predicted_power, 0.8)
predicting_power(predicted_power, 0.9)

# # Plot estimation error (median and IQR) for each scenario
# G2 <- error_summary %>%
#   mutate(group = paste0(as.character(day_plans), as.character(N_swabs_per_day)) %>% as.factor()) %>%
#   filter(N %in% c(40, 120, 200, 280, 360)) %>%
#   ggplot(aes(x = as.factor(N), y = med_error, col = day_plans, 
#              group = group, linetype = N_swabs_per_day,
#              shape = N_swabs_per_day)) +
#   geom_point(position = position_dodge(width = 0.7), size = 2.5, alpha = 0.8) +
#   geom_errorbar(aes(ymin = q1_error , ymax = q3_error ),
#                 position = position_dodge(width = 0.7), width = 0, alpha = 0.8, linewidth = 0.75) +
#   facet_wrap(trt_effect_comp ~ ., ncol = 3) +
#   theme_bw(base_size = 13) +
#   geom_hline(yintercept = 0, linetype = 'dashed', col = 'red' , linewidth = 0.7, alpha = 0.5) +
#   scale_color_manual(values = (c("#3D90D7", "#B03052")), name = "") +
#   scale_linetype_manual(values = c("dashed", "solid"), name = "") +
#   scale_shape_manual(values = c(16,17), name = "") +
#   xlab("Number of patients per arm") +
#   ylab("Estimation error") +
#   ylim(-0.25,0.25)
# G2
# ggsave(here("Plots", "simulation_estimation_error.png"), G2, width = 9, height = 8, units = 'in', dpi = 300)

library(ggpubr)
library(gridExtra)
library(gtable)

table_for_plot <- function(table){
  G <- table %>%
    tableGrob(rows = NULL, 
              theme = ttheme_minimal()
    ) %>% gtable_add_grob(grobs = grid.rect(gp = gpar(fill = NA, lwd = 2)),
                          t = 2, b = nrow(.), l = 1, r = ncol(.)) %>%
    gtable_add_grob(grobs = grid.rect(gp = gpar(fill = NA, lwd = 2)),
                    t = 1, l = 1, r = ncol(.)) %>%
    gtable_add_grob(grobs = grid.rect(gp = gpar(fill = NA, lwd = 2)),
                    t = 1, l = 1, r = 1, b = nrow(.)) %>%
    gtable_add_grob(grobs = grid.rect(gp = gpar(fill = NA, lwd = 2)),
                    t = 1, l = 2, r = 2, b = nrow(.)) %>%
    gtable_add_grob(grobs = grid.rect(gp = gpar(fill = NA, lwd = 2)),
                    t = 1, l = 4, r = 4, b = nrow(.)) %>%
    as_ggplot()
  
  G
} 

combined_plot <-  ggarrange(
  G3 %>%
    annotate_figure(top = text_grob("A) Statistical power", face = "bold", 
                                    size = 15, hjust = 2.25)),
  predicting_power(predicted_power, 0.8)  %>%
    table_for_plot() %>%
    annotate_figure(top = text_grob("B) Sample size required for 80% power", face = "bold", 
                                    size = 15, hjust = 1.125)),
  predicting_power(predicted_power, 0.9) %>%
    table_for_plot() %>%
    annotate_figure(top = text_grob("C) Sample size required for 90% power", face = "bold", 
                                    size = 15, hjust = 1.125)),
  ncol = 1,
  heights = c(3, 1,1)
)

pdf(here("Plots", "simulation_final_res.pdf"), width = 9.5, height = 9.5)
combined_plot
dev.off()

png(here("Plots", "simulation_final_res.png"), width = 9.5, height = 9.5, units = 'in', res = 350)
combined_plot
dev.off()





