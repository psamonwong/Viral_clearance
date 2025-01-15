res_neg_control_all <- res_all %>% filter(trt_control == 1) %>%
  mutate(trt_effect_comp_lab = as.factor(paste0("Treatment effect = ", (trt_effect_comp-1)*100, "%"))) %>%
  group_by(N, trt_effect_comp) %>%
  arrange(Median) %>%
  mutate(n_per_group = n(),
         sim_k2 = 1:n_per_group,
         median_effect = (Median),
         diff_estimates = median_effect - log(trt_effect_comp)) %>%
  ungroup()

res_neg_control_all$trt_effect_comp_lab <- factor(res_neg_control_all$trt_effect_comp_lab,
                                          levels = (paste0("Treatment effect = ", seq(20,100,20), "%")))


# res_neg_control_summarise <- res_neg_control_all %>%
#   group_by(N, trt_effect_comp) %>%
#   summarise(med_diff = median(diff_estimates),
#             low_diff = quantile(diff_estimates, 0.25),
#             up_diff = quantile(diff_estimates, 0.75),
#             range_diff = up_diff - low_diff,
#             lab = paste0(round(med_diff, 1), " [IQR: ",
#             round(low_diff, 1), " to ",
#             round(up_diff,1), "]"))


# ggplot(res_neg_control_all) +
#   geom_point(aes(x = diff_estimates, y = sim_k2)) +
#   facet_grid(N ~ trt_effect_comp) +
#   xlim(-2,2) +
#   geom_vline(xintercept = 0, col = 'red', linetype = 'dashed',
#              linewidth = 0.75) +
#   theme_bw() +
#   geom_text(data = res_neg_control_summarise,
#             x = 1.5, y = 90, aes(label = round(range_diff, 2)))


#Precision
# ggplot(res_neg_control_summarise,
#        aes(x = N, y = range_diff, col = as.factor(trt_effect_comp))) +
#   geom_point() +
#   geom_line() +
#   theme_bw() +
#   ylim(0,0.4)


#Accuracies
# ggplot(res_neg_control_summarise,
#        aes(x = N, y = med_diff, col = as.factor(trt_effect_comp))) +
#   geom_point() +
#   geom_line() +
#   theme_bw() +
#   ylim(-0.5, 0.5)



G <- ggplot(res_neg_control_all,
       aes(x = as.factor(N), y = diff_estimates)) +
  geom_jitter(alpha = 0.15, width = 0.2) +
  geom_boxplot(fill = '#6A80B9', alpha = 0.8, lwd = 0.5,
               outlier.shape = NA) +
  facet_wrap(. ~ trt_effect_comp_lab) +
  geom_hline(yintercept = 0, col = 'red', linewidth = 0.75, linetype = 'dashed') +
  theme_bw(base_size = 12) +
  ylab("Estimation error") +
  xlab("Number of patients per arm")
G


png("Plots/estimation_error.png", width = 8, height = 6, unit = "in", res = 350)
G
dev.off()  
  
  
