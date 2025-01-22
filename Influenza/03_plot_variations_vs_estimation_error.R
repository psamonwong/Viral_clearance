load('Rout/sim_settings_sigmasq_u12.RData')
###################################################################################################
library(ggplot2)
library(dplyr)
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
res_sigmasq_u2 <- res_all %>% filter(k_sigmasq_u_1 == 1) %>%
  mutate(trt_effect_comp_lab = as.factor(paste0("Treatment effect = ", (trt_effect_comp-1)*100, "%"))) %>%
  group_by(N, k_sigmasq_u_2, trt_effect_comp) %>%
  arrange(Median) %>%
  mutate(n_per_group = n(),
         sim_k2 = 1:n_per_group,
         median_effect = (Median),
         diff_estimates = median_effect - log(trt_effect_comp)) %>%
  ungroup()

res_sigmasq_u2$trt_effect_comp_lab <- factor(res_sigmasq_u2$trt_effect_comp_lab,
                                             levels = (paste0("Treatment effect = ", seq(20,100,20), "%")))

res_sigmasq_u2$k_sigmasq_u_2 <- as.factor(res_sigmasq_u2$k_sigmasq_u_2)
res_sigmasq_u2$facet <- c("0.5x", "1.0x", "2.0x")[as.numeric(res_sigmasq_u2$k_sigmasq_u_2)]

str(res_sigmasq_u2)

diff_sigmasq_u2 <-  res_sigmasq_u2 %>%
  group_by(N, trt_effect_comp_lab,k_sigmasq_u_2) %>%
  summarise(med = median(diff_estimates),
            low_IQR = quantile(diff_estimates, 0.25),
            up_IQR = quantile(diff_estimates, 0.75),
            min = min(diff_estimates),
            max = max(diff_estimates)) %>%
  mutate(N = as.factor(N),
         k_sigmasq_u_2 = as.factor(k_sigmasq_u_2)) %>%
  as.data.frame()

levels(diff_sigmasq_u2$k_sigmasq_u_2) <- c("0.5x", "1.0x", "2.0x")


G1 <- ggplot(diff_sigmasq_u2, aes(x = N, y = med, col = k_sigmasq_u_2)) +
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
  ggtitle(expression("C) Varied inter-individual variation on the slope (" * sigma[theta[2]]^2 * ")")) +
  scale_color_manual(values = c("#DDA853", "#DA498D", "#69247C"), name = "") +
  theme(plot.margin = unit(c(0,0.1,0,0), 'lines'))

png("Plots/estimation_error_sigmasq_u2.png", width = 8, height = 6, unit = "in", res = 350)
G1
dev.off()    

###################################################################################################
# Interindividual variation on intercept (sigmasq_u1)
res_sigmasq_u1 <- res_all %>% filter(k_sigmasq_u_2 == 1) %>%
  mutate(trt_effect_comp_lab = as.factor(paste0("Treatment effect = ", (trt_effect_comp-1)*100, "%"))) %>%
  group_by(N, k_sigmasq_u_1, trt_effect_comp) %>%
  arrange(Median) %>%
  mutate(n_per_group = n(),
         sim_k2 = 1:n_per_group,
         median_effect = (Median),
         diff_estimates = median_effect - log(trt_effect_comp)) %>%
  ungroup()

res_sigmasq_u1$trt_effect_comp_lab <- factor(res_sigmasq_u1$trt_effect_comp_lab,
                                             levels = (paste0("Treatment effect = ", seq(20,100,20), "%")))

res_sigmasq_u1$k_sigmasq_u_1 <- as.factor(res_sigmasq_u1$k_sigmasq_u_1)
res_sigmasq_u1$facet <- c("0.5x", "1.0x", "2.0x")[as.numeric(res_sigmasq_u1$k_sigmasq_u_1)]

diff_sigmasq_u1 <-  res_sigmasq_u1 %>%
  group_by(N, trt_effect_comp_lab,k_sigmasq_u_1) %>%
  summarise(med = median(diff_estimates),
            low_IQR = quantile(diff_estimates, 0.25),
            up_IQR = quantile(diff_estimates, 0.75),
            min = min(diff_estimates),
            max = max(diff_estimates)) %>%
  mutate(N = as.factor(N),
         k_sigmasq_u_1 = as.factor(k_sigmasq_u_1)) %>%
  as.data.frame()

levels(diff_sigmasq_u1$k_sigmasq_u_1) <- c("0.5x", "1.0x", "2.0x")

G2 <- ggplot(diff_sigmasq_u1, aes(x = N, y = med, col = k_sigmasq_u_1)) +
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
  ggtitle(expression("B) Varied inter-individual variation on the intercept (" * sigma[theta[1]]^2 * ")")) +
  scale_color_manual(values = c("#DDA853", "#DA498D", "#69247C"), name = "") +
  theme(plot.margin = unit(c(0,0.1,0,0), 'lines'))

G2

png("Plots/estimation_error_sigmasq_u1.png", width = 8, height = 6, unit = "in", res = 350)
G2
dev.off()  


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



res_logvl <- res_all %>% filter(trt_control == 1) %>%
  mutate(trt_effect_comp_lab = as.factor(paste0("Treatment effect = ", (trt_effect_comp-1)*100, "%"))) %>%
  group_by(N, k_sigma_logvl, trt_effect_comp) %>%
  arrange(Median) %>%
  mutate(n_per_group = n(),
         sim_k2 = 1:n_per_group,
         median_effect = (Median),
         diff_estimates = median_effect - log(trt_effect_comp)) %>%
  ungroup()

res_logvl$trt_effect_comp_lab <- factor(res_logvl$trt_effect_comp_lab,
                                             levels = (paste0("Treatment effect = ", seq(20,100,20), "%")))

res_logvl$k_sigma_logvl <- as.factor(res_logvl$k_sigma_logvl)
res_logvl$facet <- c("0.5x", "1.0x", "2.0x")[as.numeric(res_logvl$k_sigma_logvl)]


diff_sigma_logvl <-  res_logvl %>%
  group_by(N, trt_effect_comp_lab,k_sigma_logvl) %>%
  summarise(med = median(diff_estimates),
            low_IQR = quantile(diff_estimates, 0.25),
            up_IQR = quantile(diff_estimates, 0.75),
            min = min(diff_estimates),
            max = max(diff_estimates)) %>%
  mutate(N = as.factor(N),
         k_sigma_logvl = as.factor(k_sigma_logvl)) %>%
  as.data.frame()

levels(diff_sigma_logvl$k_sigma_logvl) <- c("0.5x", "1.0x", "2.0x")

G3 <- ggplot(diff_sigma_logvl, aes(x = N, y = med, col = k_sigma_logvl)) +
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
  ggtitle(expression("A) Varied observation error (" * sigma[logvl] * ")")) +
  scale_color_manual(values = c("#DDA853", "#DA498D", "#69247C"), name = "") +
  theme(plot.margin = unit(c(0,0.1,0,0), 'lines'))

G3


png("Plots/estimation_error_sigma_logvl.png", width = 8, height = 6, unit = "in", res = 350)
G3
dev.off()  




library(ggpubr)
library(grid)

plot_combined <- ggarrange(G3, G2, G1, nrow = 3, align = "v",
                           common.legend = T, legend = "right")


png("Plots/estimation_error.png", width = 11, height = 8, unit = "in", res = 350)
print(annotate_figure(plot_combined, bottom = textGrob("Number of patients per arm", vjust = 0.5, gp = gpar(cex = 1.2, fontface="bold")),
                      left = textGrob("Estimation error", rot = 90, gp = gpar(cex = 1.2, fontface="bold"))))

dev.off()  

?ggarrange
