load('Rout/sim_settings.RData')

library(ggplot2)
library(dplyr)

res <- array(dim = c(nrow(sim_settings), 3))

for(i in 1:nrow(sim_settings)){
  res_data <- read.csv(paste('sims_out/sim_out',i,'.csv',sep=''))
  res[i,] <- quantile(res_data$trt_effect, c(0.5, 0.025, 0.975))
}
colnames(res) <- c("Median", "Lower", "Upper")
res_all <- cbind(sim_settings, res)
###################################################################################################
res_neg_control <- res_all %>% filter(trt_control == 1) %>%
  mutate(sig = if_else(Lower > 0, T, F)) %>%
  group_by(N, trt_effect_comp) %>%
  summarise(total = n(),
            n = sum(sig),
            power = n/total) %>%
  mutate(trt_effect_comp = as.factor(paste0("Effect size = ", (trt_effect_comp-1)*100, "%")))

res_neg_control$trt_effect_comp <- factor(res_neg_control$trt_effect_comp,
                                          levels = rev(paste0("Effect size = ", seq(20,100,20), "%")))


A <- ggplot(res_neg_control, aes(x = N, y = power, col = trt_effect_comp)) +
  geom_point(size = 3.5, alpha = 0.8) +
  geom_line(linewidth = 0.8, linetype = "dashed") +
  theme_bw(base_size = 8) +
  ylab("Power") +
  xlab("Number of patients per arm") +
  scale_x_continuous(breaks = seq(40,240,40)) +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
  scale_color_manual(values = c("#1230AE", "#6C48C5", "#C68FE6",
                                "#D8A25E", "#A04747"),
                     name = "") +
  ggtitle("A) Baseline clearance kinetics") +
  geom_hline(yintercept = 0.8, linetype = "dashed", col = "red")
A
###################################################################################################
# res_pos_control <- res_all %>% filter(trt_control == 1.6) %>%
#   mutate(sig = if_else(Upper < 0, T, F)) %>%
#   group_by(N, trt_effect_comp) %>%
#   summarise(total = n(),
#             n = sum(sig),
#             power = n/total) %>%
#   mutate(trt_effect_comp = as.factor(paste0("Effect size = ", (trt_effect_comp-1)*100, "%")))
# 
# res_pos_control$trt_effect_comp <- factor(res_pos_control$trt_effect_comp,
#                                           levels = rev(paste0("Effect size = ", seq(20,100,20), "%")))
# 
# 
# B <- ggplot(res_pos_control, aes(x = N, y = power, col = trt_effect_comp)) +
#   geom_point(size = 2.5, alpha = 0.8) +
#   geom_line(linewidth = 0.75, linetype = "dashed") +
#   theme_bw(base_size = 8) +
#   ylab("Proportion of rejecting H0 (inferior)") +
#   xlab("Number of patients per arm") +
#   scale_x_continuous(breaks = seq(40,240,40)) +
#   scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
#   scale_color_manual(values = c("#1230AE", "#6C48C5", "#C68FE6",
#                                 "#D8A25E", "#A04747"),
#                      name = "") +
#   ggtitle("Comparing against positive control with 60% effect") +
#   geom_hline(yintercept = 0.8, linetype = "dashed", col = "red")
# B
# ###################################################################################################
# library(ggpubr)
# library(cowplot)
# 
# 
# plot_combinded <- ggarrange(A,B, ncol = 2, common.legend = T, labels = "AUTO", legend = "right",
#           align = "hv")
# 
# png("Plots/Power_analysis.png", width = 7, height = 3.5, unit = "in", res = 350)
# plot_combinded
# dev.off()
###################################################################################################
load('Rout/sim_settings_sigmasq_u12.RData')
res <- array(dim = c(nrow(sim_settings), 3))

for(i in 1:nrow(sim_settings)){
  res_data <- read.csv(paste('sims_out_sigmasq_u12/sim_out_sigmasq_u12_',i,'.csv',sep=''))
  res[i,] <- quantile(res_data$trt_effect, c(0.5, 0.025, 0.975))
}
colnames(res) <- c("Median", "Lower", "Upper")
res_all <- cbind(sim_settings, res)
res_1 <- res_all
###################################################################################################
load('Rout/sim_settings_sigmasq_u12_supplement.RData')
res <- array(dim = c(nrow(sim_settings), 3))

for(i in 1:nrow(sim_settings)){
  res_data <- read.csv(paste('sims_out_sigmasq_u12_supplement/sim_out_sigmasq_u12_supplement',i,'.csv',sep=''))
  res[i,] <- quantile(res_data$trt_effect, c(0.5, 0.025, 0.975))
}
colnames(res) <- c("Median", "Lower", "Upper")
res_all <- cbind(sim_settings, res)
res_2 <- res_all
###################################################################################################
res_all <- rbind(res_1, res_2)
###################################################################################################
# Interindividual variation on slope (sigmasq_u2)
res_sigmasq_u2 <- res_all %>%
  filter(k_sigmasq_u_1 == 1) %>%
  mutate(sig = if_else(Lower > 0, T, F)) %>%
  group_by(N, k_sigmasq_u_2, trt_effect_comp) %>%
  summarise(total = n(),
            n = sum(sig),
            power = n/total) %>%
  mutate(trt_effect_comp = as.factor(paste0("Effect size = ", (trt_effect_comp-1)*100, "%")))

res_sigmasq_u2$trt_effect_comp <- factor(res_sigmasq_u2$trt_effect_comp,
                                         levels = rev(paste0("Effect size = ", seq(20,100,20), "%")))

res_sigmasq_u2$k_sigmasq_u_2 <- as.factor(res_sigmasq_u2$k_sigmasq_u_2)
res_sigmasq_u2$facet <- c("0.5x", "1.0x", "2.0x")[as.numeric(res_sigmasq_u2$k_sigmasq_u_2)]

res_sigmasq_u2 <- res_sigmasq_u2 %>% filter(facet != "1.0x")

G1 <- ggplot(res_sigmasq_u2, aes(x = N, y = power, col = trt_effect_comp)) +
  geom_point(size = 2, alpha = 0.8) +
  facet_wrap(.~facet) +
  geom_line(linewidth = 0.5, linetype = "dashed") +
  theme_bw(base_size = 10) +
  ylab("") +
  xlab("") +
  scale_x_continuous(breaks = seq(40,240,40)) +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
  scale_color_manual(values = c("#1230AE", "#6C48C5", "#C68FE6",
                                "#D8A25E", "#A04747"),
                     name = "") +
  ggtitle(expression("C) Varied inter-individual variation on the slope (" * sigma[theta[2]]^2 * ")")) +
  geom_hline(yintercept = 0.8, linetype = "dashed", col = "red") +
  theme(plot.margin = unit(c(0,0.2,0,0), 'lines'))
G1
###################################################################################################
# Interindividual variation on intercept (sigmasq_u1)
res_sigmasq_u1 <- res_all %>%
  filter(k_sigmasq_u_2 == 1) %>%
  mutate(sig = if_else(Lower > 0, T, F)) %>%
  group_by(N, k_sigmasq_u_1, trt_effect_comp) %>%
  summarise(total = n(),
            n = sum(sig),
            power = n/total) %>%
  mutate(trt_effect_comp = as.factor(paste0("Effect size = ", (trt_effect_comp-1)*100, "%")))

res_sigmasq_u1$trt_effect_comp <- factor(res_sigmasq_u1$trt_effect_comp,
                                         levels = rev(paste0("Effect size = ", seq(20,100,20), "%")))

res_sigmasq_u1$k_sigmasq_u_1 <- as.factor(res_sigmasq_u1$k_sigmasq_u_1)
res_sigmasq_u1$facet <- c("0.5x", "1.0x", "2.0x")[as.numeric(res_sigmasq_u1$k_sigmasq_u_1)]

res_sigmasq_u1 <- res_sigmasq_u1 %>% filter(facet != "1.0x")


G2 <- ggplot(res_sigmasq_u1, aes(x = N, y = power, col = trt_effect_comp)) +
  geom_point(size = 2, alpha = 0.8) +
  facet_wrap(.~facet) +
  geom_line(linewidth = 0.5, linetype = "dashed") +
  theme_bw(base_size = 10) +
  ylab("") +
  xlab("") +
  scale_x_continuous(breaks = seq(40,240,40)) +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
  scale_color_manual(values = c("#1230AE", "#6C48C5", "#C68FE6",
                                "#D8A25E", "#A04747"),
                     name = "") +
  ggtitle(expression("B) Varied inter-individual variation on the intercept (" * sigma[theta[1]]^2 * ")")) +
  geom_hline(yintercept = 0.8, linetype = "dashed", col = "red") +
  theme(plot.margin = unit(c(0,0.2,0,0), 'lines'))

G2
###################################################################################################
###################################################################################################
load('Rout/sim_settings_sigma_logvl.RData')

res <- array(dim = c(nrow(sim_settings), 3))

for(i in 1:nrow(sim_settings)){
  res_data <- read.csv(paste('sims_out_sigma_logvl/sim_out_sigma_logvl_',i,'.csv',sep=''))
  res[i,] <- quantile(res_data$trt_effect, c(0.5, 0.025, 0.975))
}
colnames(res) <- c("Median", "Lower", "Upper")
res_all <- cbind(sim_settings, res)
res_1 <- res_all

###################################################################################################
load('Rout/sim_settings_sigma_logvl_supplement.RData')

res <- array(dim = c(nrow(sim_settings), 3))

for(i in 1:nrow(sim_settings)){
  res_data <- read.csv(paste('sims_out_sigma_logvl_supplement/sim_out_sigma_logvl_supplement_',i,'.csv',sep=''))
  res[i,] <- quantile(res_data$trt_effect, c(0.5, 0.025, 0.975))
}
colnames(res) <- c("Median", "Lower", "Upper")
res_all <- cbind(sim_settings, res)
res_2 <- res_all

res_all <- rbind(res_1, res_2)
###################################################################################################
res_logvl <- res_all %>% filter(trt_control == 1) %>%
  mutate(sig = if_else(Lower > 0, T, F)) %>%
  group_by(N, k_sigma_logvl, trt_effect_comp) %>%
  summarise(total = n(),
            n = sum(sig),
            power = n/total) %>%
  mutate(trt_effect_comp = as.factor(paste0("Effect size = ", (trt_effect_comp-1)*100, "%")))

res_logvl$trt_effect_comp <- factor(res_logvl$trt_effect_comp,
                                    levels = rev(paste0("Effect size = ", seq(20,100,20), "%")))

res_logvl$k_sigma_logvl <- as.factor(res_logvl$k_sigma_logvl)
res_logvl$facet <- c("0.5x", "1.0x", "2.0x")[as.numeric(res_logvl$k_sigma_logvl)]
res_logvl <- res_logvl %>% filter(facet != "1.0x")


G3 <- ggplot(res_logvl, aes(x = N, y = power, col = trt_effect_comp)) +
  geom_point(size = 2, alpha = 0.8) +
  facet_wrap(.~as.factor(facet)) +
  geom_line(linewidth = 0.5, linetype = "dashed") +
  theme_bw(base_size = 10) +
  ylab("") +
  xlab("") +
  scale_x_continuous(breaks = seq(40,240,40)) +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
  scale_color_manual(values = c("#1230AE", "#6C48C5", "#C68FE6",
                                "#D8A25E", "#A04747"),
                     name = "") +
  ggtitle(expression("D) Varied observation error (" * sigma[VL]^2 * ")")) +
  geom_hline(yintercept = 0.8, linetype = "dashed", col = "red") +
  theme(plot.margin = unit(c(0,0.2,0,0), 'lines'))
G3
###################################################################################################
library(ggpubr)
library(grid)

combined_plot <- ggarrange(G2, G1, G3, 
                           common.legend = T, ncol = 1,
                           legend = "right",
                           align = "hv")


combined_plot <- annotate_figure( combined_plot, 
                 bottom = textGrob("Number of patients per arm", vjust = 0.5, gp = gpar(cex = 1.2, fontface="bold")),
                 left = textGrob("Power", rot = 90, gp = gpar(cex = 1.2, fontface="bold")))


combined_plot_all <- ggarrange(A + theme_bw(base_size = 14) +
            theme(legend.position = "none",
                  axis.title = element_text(face = "bold")) , 
          combined_plot)
combined_plot_all

plot_name <- paste0('Plots/variations_vs_power.png')

png(plot_name, width = 11,height = 6.5, units = 'in', res = 350)
print(combined_plot_all)
dev.off()





