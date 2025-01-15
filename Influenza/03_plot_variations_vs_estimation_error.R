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

G1 <- ggplot(res_sigmasq_u2,
       aes(x = as.factor(N), y = diff_estimates)) +
  geom_jitter(alpha = 0.15, width = 0.2) +
  geom_boxplot(fill = '#6A80B9', alpha = 0.8, lwd = 0.5,
               outlier.shape = NA) +
  facet_grid(trt_effect_comp_lab ~ facet) +
  geom_hline(yintercept = 0, col = 'red', linewidth = 0.75, linetype = 'dashed') +
  theme_bw(base_size = 12) +
  ylab("Estimation error") +
  xlab("Number of patients per arm") +
  ggtitle(expression("Inter-individual variation on the slope (" * sigma[u2]^2 * ")")) +
  ylim(-1,1)
  
G1

png("Plots/estimation_error_sigmasq_u2.png", width = 8, height = 8, unit = "in", res = 350)
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

G2 <- ggplot(res_sigmasq_u1,
             aes(x = as.factor(N), y = diff_estimates)) +
  geom_jitter(alpha = 0.15, width = 0.2) +
  geom_boxplot(fill = '#6A80B9', alpha = 0.8, lwd = 0.5,
               outlier.shape = NA) +
  facet_grid(trt_effect_comp_lab ~ facet) +
  geom_hline(yintercept = 0, col = 'red', linewidth = 0.75, linetype = 'dashed') +
  theme_bw(base_size = 12) +
  ylab("Estimation error") +
  xlab("Number of patients per arm") +
  ggtitle(expression("Inter-individual variation on the intercept (" * sigma[u1]^2 * ")")) +
  ylim(-1,1)
G2

png("Plots/estimation_error_sigmasq_u1.png", width = 8, height = 8, unit = "in", res = 350)
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

G3 <- ggplot(res_logvl,
             aes(x = as.factor(N), y = diff_estimates)) +
  geom_jitter(alpha = 0.15, width = 0.2) +
  geom_boxplot(fill = '#6A80B9', alpha = 0.8, lwd = 0.5,
               outlier.shape = NA) +
  facet_grid(trt_effect_comp_lab ~ facet) +
  geom_hline(yintercept = 0, col = 'red', linewidth = 0.75, linetype = 'dashed') +
  theme_bw(base_size = 12) +
  ylab("Estimation error") +
  xlab("Number of patients per arm") +
  ggtitle(expression("Swabbing variation (" * sigma[logvl] * ")")) +
  ylim(-1,1)
G3

png("Plots/estimation_error_sigma_logvl.png", width = 8, height = 8, unit = "in", res = 350)
G3
dev.off()  
