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
# Sampling day
res_days <- res_all %>%
  filter(N_swabs_per_day == 2) %>%
  mutate(sig = if_else(Lower > 0, T, F)) %>%
  group_by(N, day_plans, trt_effect_comp) %>%
  summarise(total = n(),
            n = sum(sig),
            power = n/total) %>%
  mutate(trt_effect_comp = as.factor(paste0("Effect size = ", (trt_effect_comp-1)*100, "%")))

res_days$trt_effect_comp <- factor(res_days$trt_effect_comp,
                                          levels = rev(paste0("Effect size = ", seq(20,100,20), "%")))

res_days$lab <- NA
res_days$lab[res_days$day_plans == "0,1,2,3,4,5"] <- "Everyday \n(n=6)"
res_days$lab[res_days$day_plans == "0,2,5"] <- "Everyday other day \n(n=3)"
res_days$lab[res_days$day_plans == "0,5"] <- "First and last day \n(n=2)"


#res_sigmasq_u2$k_sigmasq_u_2 <- as.factor(res_sigmasq_u2$k_sigmasq_u_2)
#res_sigmasq_u2$facet <- c("0.5x", "1.0x", "2.0x")[as.numeric(res_sigmasq_u2$k_sigmasq_u_2)]

G1 <- ggplot(res_days, aes(x = N, y = power, col = trt_effect_comp)) +
  geom_point(size = 2, alpha = 0.8) +
  facet_wrap(.~lab) +
  geom_line(linewidth = 0.5, linetype = "dashed") +
  theme_bw(base_size = 13) +
  ylab("Power") +
  xlab("Number of patients per arm") +
  scale_x_continuous(breaks = seq(40,240,40)) +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
  scale_color_manual(values = c("#1230AE", 
                                "#6C48C5", 
                                "#C68FE6",
                                "#D8A25E", 
                                "#A04747"
                                ),
                     name = "") +
  geom_hline(yintercept = 0.8, linetype = "dashed", col = "red")  +
  ggtitle("A) Varied sampling schedule") +
  theme(axis.title = element_text(face = "bold"))
G1
###################################################################################################
# Swabs day
res_swabs <- res_all %>%
  filter(day_plans == "0,1,2,3,4,5") %>%
  mutate(sig = if_else(Lower > 0, T, F)) %>%
  group_by(N, N_swabs_per_day, trt_effect_comp) %>%
  summarise(total = n(),
            n = sum(sig),
            power = n/total) %>%
  mutate(trt_effect_comp = as.factor(paste0("Effect size = ", (trt_effect_comp-1)*100, "%")),
         lab = if_else(N_swabs_per_day == 1, paste0(N_swabs_per_day, " swab/day"), paste0(N_swabs_per_day, " swabs/day")))

res_swabs$trt_effect_comp <- factor(res_swabs$trt_effect_comp,
                                   levels = rev(paste0("Effect size = ", seq(20,100,20), "%")))

#res_days$lab[res_days$day_plans == "0,1,2,3,4,5"] <- "Everyday \n(n=6)"
#res_days$lab[res_days$day_plans == "0,2,5"] <- "Everyday other day \n(n=3)"
#res_days$lab[res_days$day_plans == "0,5"] <- "Fist and last day \n(n=2)"


#res_sigmasq_u2$k_sigmasq_u_2 <- as.factor(res_sigmasq_u2$k_sigmasq_u_2)
#res_sigmasq_u2$facet <- c("0.5x", "1.0x", "2.0x")[as.numeric(res_sigmasq_u2$k_sigmasq_u_2)]

G2 <- ggplot(res_swabs, aes(x = N, y = power, col = trt_effect_comp)) +
  geom_point(size = 2, alpha = 0.8) +
  facet_wrap(.~lab) +
  geom_line(linewidth = 0.5, linetype = "dashed") +
  theme_bw(base_size = 13) +
  ylab("Power") +
  xlab("Number of patients per arm") +
  scale_x_continuous(breaks = seq(40,240,40)) +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
  scale_color_manual(values = c("#1230AE", 
                                "#6C48C5", 
                                "#C68FE6",
                                "#D8A25E", 
                                "#A04747"
  ),
  name = "") +
  geom_hline(yintercept = 0.8, linetype = "dashed", col = "red") +
  ggtitle("B) Varied number of swabs per day") +
  theme(axis.title = element_text(face = "bold"))
G2
####################################################################
library(ggpubr)
library(grid)

plot_combined <- ggarrange(G1, G2, nrow = 2, align = "v",
                           common.legend = T, legend = "right")


png("Plots/sampling_power.png", width = 10, height = 8, unit = "in", res = 350)
plot_combined
dev.off()  





























# Sampling frequency
res_sampling_fq <- res_all %>%
  filter(day_plans == c("0,1,2,3,4,5")) %>%
  mutate(sig = if_else(Lower > 0, T, F)) %>%
  group_by(N, N_swabs_per_day, trt_effect_comp) %>%
  summarise(total = n(),
            n = sum(sig),
            power = n/total) %>%
  mutate(trt_effect_comp = as.factor(paste0("Effect size = ", (trt_effect_comp-1)*100, "%")))

res_sampling_fq$trt_effect_comp <- factor(res_sampling_fq$trt_effect_comp,
                                         levels = rev(paste0("Effect size = ", seq(20,100,20), "%")))

#res_sigmasq_u1$k_sigmasq_u_1 <- as.factor(res_sigmasq_u1$k_sigmasq_u_1)
#res_sigmasq_u1$facet <- c("0.5x", "1.0x", "2.0x")[as.numeric(res_sigmasq_u1$k_sigmasq_u_1)]

 ggplot(res_sampling_fq, aes(x = N, y = power, col = trt_effect_comp)) +
  geom_point(size = 2, alpha = 0.8) +
  facet_wrap(.~N_swabs_per_day) +
  geom_line(linewidth = 0.5, linetype = "dashed") +
  theme_bw(base_size = 10) +
  ylab("") +
  xlab("") +
  scale_x_continuous(breaks = seq(40,240,40)) +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
   scale_color_manual(values = c("#1230AE", 
                                 #"#6C48C5", 
                                 "#C68FE6",
                                 #"#D8A25E", 
                                 "#A04747"
   ),
   name = "") +
  #ggtitle(expression("A) Inter-individual variation on the intercept (" * sigma[theta[1]]^2 * ")")) +
  geom_hline(yintercept = 0.8, linetype = "dashed", col = "red")

G2

res_N_swabs_etimation_error <- res_all %>% 
  filter(day_plans == c("0,1,2,3,4,5")) %>%
  mutate(trt_effect_comp_lab = as.factor(paste0("Effect size = ", (trt_effect_comp-1)*100, "%"))) %>%
  group_by(N, N_swabs_per_day, trt_effect_comp) %>%
  arrange(Median) %>%
  mutate(n_per_group = n(),
         sim_k2 = 1:n_per_group,
         median_effect = (Median),
         diff_estimates = median_effect - log(trt_effect_comp)) %>%
  ungroup()

res_N_swabs_etimation_error$trt_effect_comp_lab <- factor(res_N_swabs_etimation_error$trt_effect_comp_lab,
                                                        levels = (paste0("Effect size = ", seq(20,100,20), "%")))

#res_sigmasq_u2$k_sigmasq_u_2 <- as.factor(res_sigmasq_u2$k_sigmasq_u_2)
#res_sigmasq_u2$facet <- c("0.5x", "1.0x", "2.0x")[as.numeric(res_sigmasq_u2$k_sigmasq_u_2)]

#str(res_sigmasq_u2)

diff_N_swabs_estimation_error <-  res_N_swabs_etimation_error %>%
  group_by(N, trt_effect_comp_lab,N_swabs_per_day) %>%
  summarise(med = median(diff_estimates),
            low_IQR = quantile(diff_estimates, 0.25),
            up_IQR = quantile(diff_estimates, 0.75),
            min = min(diff_estimates),
            max = max(diff_estimates)) %>%
  mutate(N = as.factor(N),
         N_swabs_per_day = as.factor(N_swabs_per_day)) %>%
  as.data.frame()

#levels(diff_sigmasq_u2$k_sigmasq_u_2) <- c("0.5x", "1.0x", "2.0x")


ggplot(diff_N_swabs_estimation_error, aes(x = N, y = med, col = N_swabs_per_day)) +
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
  #ggtitle(expression("C) Varied inter-individual variation on the slope (" * sigma[theta[2]]^2 * ")")) +
  scale_color_manual(values = c("#DDA853", "#DA498D", "#69247C"), name = "") +
  theme(plot.margin = unit(c(0,0.1,0,0), 'lines'))













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
                       left = textGrob("Power", rot = 90, gp = gpar(cex = 1.2, fontface="bold"))))
dev.off()
