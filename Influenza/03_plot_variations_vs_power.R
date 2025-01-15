load('Rout/sim_settings_sigmasq_u12.RData')
###################################################################################################
library(ggplot2)
library(dplyr)
###################################################################################################
res <- array(dim = c(nrow(sim_settings), 3))

for(i in 1:nrow(sim_settings)){
  res_data <- read.csv(paste('sims_out_sigmasq_u12/sim_out_sigmasq_u12_',i,'.csv',sep=''))
  res[i,] <- quantile(res_data$trt_effect, c(0.5, 0.025, 0.975))
}
colnames(res) <- c("Median", "Lower", "Upper")
res_all <- cbind(sim_settings, res)
###################################################################################################
# Interindividual variation on slope (sigmasq_u2)
res_sigmasq_u2 <- res_all %>%
  filter(k_sigmasq_u_1 == 1) %>%
  mutate(sig = if_else(Lower > 0, T, F)) %>%
  group_by(N, k_sigmasq_u_2, trt_effect_comp) %>%
  summarise(total = n(),
            n = sum(sig),
            power = n/total) %>%
  mutate(trt_effect_comp = as.factor(paste0("Treatment effect = ", (trt_effect_comp-1)*100, "%")))

res_sigmasq_u2$trt_effect_comp <- factor(res_sigmasq_u2$trt_effect_comp,
                                          levels = rev(paste0("Treatment effect = ", seq(20,100,40), "%")))

res_sigmasq_u2$k_sigmasq_u_2 <- as.factor(res_sigmasq_u2$k_sigmasq_u_2)
res_sigmasq_u2$facet <- c("0.5x", "1.0x", "2.0x")[as.numeric(res_sigmasq_u2$k_sigmasq_u_2)]

G1 <- ggplot(res_sigmasq_u2, aes(x = N, y = power, col = trt_effect_comp)) +
  geom_point(size = 2.75, alpha = 0.8) +
  facet_wrap(.~facet) +
  geom_line(linewidth = 0.75) +
  theme_bw(base_size = 10) +
  ylab("") +
  xlab("") +
  scale_x_continuous(breaks = seq(40,240,40)) +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
  scale_color_manual(values = c("#1230AE", "#C68FE6",
                                "#A04747"),
                     name = "") +
  ggtitle(expression("B) Inter-individual variation on the slope (" * sigma[u2]^2 * ")")) +
  geom_hline(yintercept = 0.8, linetype = "dashed", col = "red")
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
  mutate(trt_effect_comp = as.factor(paste0("Treatment effect = ", (trt_effect_comp-1)*100, "%")))

res_sigmasq_u1$trt_effect_comp <- factor(res_sigmasq_u1$trt_effect_comp,
                                         levels = rev(paste0("Treatment effect = ", seq(20,100,40), "%")))

res_sigmasq_u1$k_sigmasq_u_1 <- as.factor(res_sigmasq_u1$k_sigmasq_u_1)
res_sigmasq_u1$facet <- c("0.5x", "1.0x", "2.0x")[as.numeric(res_sigmasq_u1$k_sigmasq_u_1)]

G2 <- ggplot(res_sigmasq_u1, aes(x = N, y = power, col = trt_effect_comp)) +
  geom_point(size = 2.75, alpha = 0.8) +
  facet_wrap(.~facet) +
  geom_line(linewidth = 0.75) +
  theme_bw(base_size = 10) +
  ylab("") +
  xlab("") +
  scale_x_continuous(breaks = seq(40,240,40)) +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
  scale_color_manual(values = c("#1230AE", "#C68FE6",
                                "#A04747"),
                     name = "") +
  ggtitle(expression("A) Inter-individual variation on the intercept (" * sigma[u1]^2 * ")")) +
  geom_hline(yintercept = 0.8, linetype = "dashed", col = "red")

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
###################################################################################################
res_logvl <- res_all %>% filter(trt_control == 1) %>%
  mutate(sig = if_else(Lower > 0, T, F)) %>%
  group_by(N, k_sigma_logvl, trt_effect_comp) %>%
  summarise(total = n(),
            n = sum(sig),
            power = n/total) %>%
  mutate(trt_effect_comp = as.factor(paste0("Treatment effect = ", (trt_effect_comp-1)*100, "%")))

res_logvl$trt_effect_comp <- factor(res_logvl$trt_effect_comp,
                                          levels = rev(paste0("Treatment effect = ", seq(20,100,40), "%")))

res_logvl$k_sigma_logvl <- as.factor(res_logvl$k_sigma_logvl)
res_logvl$facet <- c("0.5x", "1.0x", "2.0x")[as.numeric(res_logvl$k_sigma_logvl)]

G3 <- ggplot(res_logvl, aes(x = N, y = power, col = trt_effect_comp)) +
  geom_point(size = 2.75, alpha = 0.8) +
  facet_wrap(.~as.factor(facet)) +
  geom_line(linewidth = 0.75) +
  theme_bw(base_size = 10) +
  ylab("") +
  xlab("") +
  scale_x_continuous(breaks = seq(40,240,40)) +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
  scale_color_manual(values = c("#1230AE", "#C68FE6",
                                "#A04747"),
                     name = "") +
  ggtitle(expression("C) Swabbing variation (" * sigma[logvl] * ")")) +
  geom_hline(yintercept = 0.8, linetype = "dashed", col = "red")
G3
###################################################################################################
library(ggpubr)
library(grid)

combined_plot <- ggarrange(G2, G1, G3, 
          common.legend = T, ncol = 1,
          legend = "right",
          align = "hv")


plot_name <- paste0('Plots/variations_vs_power.png')

png(plot_name, width = 7.5,height = 7.5, units = 'in', res = 350)
print(annotate_figure( combined_plot, 
                       bottom = textGrob("Number of patients per arm", vjust = 0.5, gp = gpar(cex = 1.2, fontface="bold")),
                       left = textGrob("Proportion of rejecting H0", rot = 90, gp = gpar(cex = 1.2, fontface="bold"))))
dev.off()
