load('Rout/sim_settings.RData')

library(ggplot2)
library(dplyr)
#library(fitdistrplus)

res <- array(dim = c(nrow(sim_settings), 3))

for(i in 1:nrow(sim_settings)){
  res_data <- read.csv(paste('sims_out/sim_out',i,'.csv',sep=''))
  res[i,] <- quantile(res_data$trt_effect, c(0.5, 0.025, 0.975))
}
colnames(res) <- c("Median", "Lower", "Upper")
res_all <- cbind(sim_settings, res)

#res_all$conditions <- as.numeric(as.factor(paste(res_all$N, res_all$day_plans, res_all$N_swabs_per_day, res_all$trt_effects, sep = "_")))
#res_all$conditions <- factor(res_all$conditions, levels = unique(res_all$conditions))
#levels(res_all$conditions) <- as.character(1:length(unique(res_all$conditions)))
#res_all$conditions <- as.numeric(res_all$conditions)

#res_all$N_swabs_per_day <- as.factor(res_all$N_swabs_per_day)
###################################################################################################
res_neg_control <- res_all %>% filter(trt_control == 1) %>%
  mutate(sig = if_else(Lower > 0, T, F)) %>%
  group_by(N, trt_effect_comp) %>%
  summarise(total = n(),
            n = sum(sig),
            power = n/total) %>%
  mutate(trt_effect_comp = as.factor(paste0("Treatment effect = ", (trt_effect_comp-1)*100, "%")))

res_neg_control$trt_effect_comp <- factor(res_neg_control$trt_effect_comp,
                                          levels = paste0("Treatment effect = ", seq(20,100,20), "%"))


A <- ggplot(res_neg_control, aes(x = N, y = power, col = trt_effect_comp)) +
  geom_point(size = 2.75, alpha = 0.8) +
  geom_line(linewidth = 0.75, linetype = "dashed") +
  theme_bw(base_size = 10) +
  ylab("Proportion of rejecting H0") +
  xlab("Number of patients per arm") +
  scale_x_continuous(breaks = seq(40,240,40)) +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
  scale_color_manual(values = c("#1230AE", "#6C48C5", "#C68FE6",
                                "#D8A25E", "#A04747"),
                     name = "") +
  ggtitle("Comparing against negative control") +
  geom_hline(yintercept = 0.8, linetype = "dashed", col = "red")
A
###################################################################################################
res_pos_control <- res_all %>% filter(trt_control == 1.6) %>%
  mutate(sig = if_else(Upper < 0, T, F)) %>%
  group_by(N, trt_effect_comp) %>%
  summarise(total = n(),
            n = sum(sig),
            power = n/total) %>%
  mutate(trt_effect_comp = as.factor(paste0("Treatment effect = ", (trt_effect_comp-1)*100, "%")))

res_pos_control$trt_effect_comp <- factor(res_pos_control$trt_effect_comp,
                                          levels = paste0("Treatment effect = ", seq(20,100,20), "%"))


B <- ggplot(res_pos_control, aes(x = N, y = power, col = trt_effect_comp)) +
  geom_point(size = 2.75, alpha = 0.8) +
  geom_line(linewidth = 0.75, linetype = "dashed") +
  theme_bw(base_size = 10) +
  ylab("Proportion of rejecting H0 (inferior)") +
  xlab("Number of patients per arm") +
  scale_x_continuous(breaks = seq(40,240,40)) +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
  scale_color_manual(values = c("#1230AE", "#6C48C5", "#C68FE6",
                                "#D8A25E", "#A04747"),
                     name = "") +
  ggtitle("Comparing against positive control with 60% effect") +
  geom_hline(yintercept = 0.8, linetype = "dashed", col = "red")
B
###################################################################################################
library(ggpubr)
library(cowplot)


plot_combinded <- ggarrange(A,B, ncol = 2, common.legend = T, labels = "AUTO", legend = "right",
          align = "hv")

png("Plots/Power_analysis.png", width = 10, height = 5, unit = "in", res = 350)
plot_combinded
dev.off()





