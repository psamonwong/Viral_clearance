###################################################################################################
library(ggplot2)
library(dplyr)
###################################################################################################
load('Rout/sim_settings_sampling_schedule.RData')
res <- array(dim = c(nrow(sim_settings), 3))

for(i in 1:nrow(sim_settings)){
  res_data <- read.csv(paste('sims_out_sampling_schedule/sim_out_sampling_schedule_',i,'.csv',sep=''))
  res[i,] <- quantile(res_data$trt_effect, c(0.5, 0.025, 0.975))
}
colnames(res) <- c("Median", "Lower", "Upper")
res_all <- cbind(sim_settings, res)

res1 <- res_all
###################################################################################################
load('Rout/sim_settings_sampling_schedule_supplement.RData')
res <- array(dim = c(nrow(sim_settings), 3))

for(i in 1:nrow(sim_settings)){
  res_data <- read.csv(paste('sims_out_sampling_schedule_supplement/sim_out_sampling_schedule_supplement',i,'.csv',sep=''))
  res[i,] <- quantile(res_data$trt_effect, c(0.5, 0.025, 0.975))
}
colnames(res) <- c("Median", "Lower", "Upper")
res_all <- cbind(sim_settings, res)
###################################################################################################
res_all <- rbind(res1, res_all)
###################################################################################################
# Sampling schedule
res_days <- res_all %>%
  filter(N_swabs_per_day == 2) %>%
  mutate(trt_effect_comp_lab = as.factor(paste0("Effect size = ", (trt_effect_comp-1)*100, "%"))) %>%
  group_by(N, day_plans, trt_effect_comp) %>%
  arrange(Median) %>%
  mutate(n_per_group = n(),
         sim_k2 = 1:n_per_group[1],
         median_effect = (Median),
         diff_estimates = median_effect - log(trt_effect_comp)) %>%
  ungroup()

res_days$trt_effect_comp_lab <- factor(res_days$trt_effect_comp_lab,
                                             levels = (paste0("Effect size = ", seq(20,100,20), "%")))

res_days$lab <- NA
res_days$lab[res_days$day_plans == "0,1,2,3,4,5"] <- "Everyday \n(n=6)"
res_days$lab[res_days$day_plans == "0,2,5"] <- "Everyday other day \n(n=3)"
res_days$lab[res_days$day_plans == "0,5"] <- "First and last day \n(n=2)"

str(res_days)

diff_days <-  res_days %>%
  group_by(N, trt_effect_comp_lab,lab) %>%
  summarise(med = median(diff_estimates),
            low_IQR = quantile(diff_estimates, 0.25),
            up_IQR = quantile(diff_estimates, 0.75),
            min = min(diff_estimates),
            max = max(diff_estimates)) %>%
  mutate(N = as.factor(N),
         lab = as.factor(lab)) %>%
  as.data.frame()

G1 <- ggplot(diff_days, aes(x = N, y = med, col = lab)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.5, alpha = 0.75) +
  geom_errorbar(aes(ymin = low_IQR, ymax = up_IQR), width = 0, linewidth = 0.6,
                position = position_dodge(width = 0.5)) +
  facet_wrap(.~trt_effect_comp_lab, nrow = 1) +
  ylim(-0.4,0.4) +
  theme_bw(base_size = 12) +
  geom_hline(yintercept = 0, col = 'red', linetype = 'dashed',
             linewidth = 0.75) +
  xlab("") +
  ylab("") +
  ggtitle(expression("A) Varied sampling schedule")) +
  scale_color_manual(values = c("#DDA853", "#DA498D", "#69247C"), name = "") +
  theme(plot.margin = unit(c(0,0.1,0,0), 'lines'))

G1

png("Plots/estimation_error_days.png", width = 8, height = 6, unit = "in", res = 350)
G1
dev.off()   
###################################################################################################
# Swabs day
res_swabs <- res_all %>%
  filter(day_plans == "0,1,2,3,4,5") %>%
  mutate(trt_effect_comp_lab = as.factor(paste0("Effect size = ", (trt_effect_comp-1)*100, "%"))) %>%
  group_by(N, N_swabs_per_day, trt_effect_comp) %>%
  arrange(Median) %>%
  mutate(n_per_group = n(),
         sim_k2 = 1:n_per_group[1],
         median_effect = (Median),
         diff_estimates = median_effect - log(trt_effect_comp)) %>%
  ungroup()
  
res_swabs$trt_effect_comp_lab <- factor(res_swabs$trt_effect_comp_lab,
                                       levels = (paste0("Effect size = ", seq(20,100,20), "%")))
  
res_swabs$lab = ifelse(res_swabs$N_swabs_per_day == 1, paste0(res_swabs$N_swabs_per_day, " swab/day"), paste0(res_swabs$N_swabs_per_day, " swabs/day"))

diff_swabs <-  res_swabs %>%
  group_by(N, trt_effect_comp_lab,lab) %>%
  summarise(med = median(diff_estimates),
            low_IQR = quantile(diff_estimates, 0.25),
            up_IQR = quantile(diff_estimates, 0.75),
            min = min(diff_estimates),
            max = max(diff_estimates)) %>%
  mutate(N = as.factor(N),
         lab = as.factor(lab)) %>%
  as.data.frame()

G2 <- ggplot(diff_swabs, aes(x = N, y = med, col = lab)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.5, alpha = 0.75) +
  geom_errorbar(aes(ymin = low_IQR, ymax = up_IQR), width = 0, linewidth = 0.6,
                position = position_dodge(width = 0.5)) +
  facet_wrap(.~trt_effect_comp_lab, nrow = 1) +
  ylim(-0.4,0.4) +
  theme_bw(base_size = 12) +
  geom_hline(yintercept = 0, col = 'red', linetype = 'dashed',
             linewidth = 0.75) +
  xlab("") +
  ylab("") +
  ggtitle(expression("B) Varied swab number")) +
  scale_color_manual(values = c("#DDA853", "#DA498D", "#69247C"), name = "") +
  theme(plot.margin = unit(c(0,0.1,0,0), 'lines'))

G2



library(ggpubr)
library(grid)

plot_combined <- ggarrange(G1, G2, nrow = 2, align = "v",
                           common.legend = F, legend = "right")


png("Plots/estimation_error_sampling.png", width = 11, height = 8, unit = "in", res = 350)
print(annotate_figure(plot_combined, bottom = textGrob("Number of patients per arm", vjust = 0.5, gp = gpar(cex = 1.2, fontface="bold")),
                      left = textGrob("Estimation error", rot = 90, gp = gpar(cex = 1.2, fontface="bold"))))

dev.off()  
