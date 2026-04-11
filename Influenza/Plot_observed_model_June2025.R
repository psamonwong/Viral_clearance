library(here)
library(tidyverse)
library(ggpubr)


model_settings <- read.csv(here("Outputs", "model_settings.csv"))

posterior_check_plot <- function(model_settings, Model,  Dmax, Type){
  ind_exp <- which(model_settings$Dmax == Dmax & model_settings$Type == Type & grepl(Model, model_settings$mod))
  load(paste0(here("Rout", paste0("model_settings_job_", ind_exp, ".RData"))))
  load(paste0(here("Rout", paste0("model_fits_job_", ind_exp, ".RData"))))
  
  preds <-  rstan::extract(out, pars = 'preds')
  adastra_dat_analysis$preds <- apply(preds$preds, 2, median) 
  adastra_dat_analysis <- adastra_dat_analysis  %>% mutate(mod = Model)
  
  
  preds2 <- cbind("ID" = adastra_dat_analysis$ID,
                  "Timepoint_ID" = adastra_dat_analysis$Timepoint_ID, 
                  as.data.frame(t(preds$preds))) 
  
  preds3 <- preds2 %>%
    group_by(ID, Timepoint_ID) %>%
    summarise(across(V1:V1000, median)) %>%
    pivot_longer(
      col = -(1:2),
      names_to = "Iteration",
      values_to = "values") %>%
    group_by(Timepoint_ID, Iteration) %>%
    summarise(Q1 = quantile(values, 0.25),
              Q2 = quantile(values, 0.5),
              Q3 = quantile(values, 0.75)
    ) %>%
    group_by(Timepoint_ID) %>%
    summarise(Q1_low = quantile(Q1, 0.025),
              Q1_med = quantile(Q1, 0.5),
              Q1_up = quantile(Q1, 0.975),
              Q2_low = quantile(Q2, 0.025),
              Q2_med = quantile(Q2, 0.5),
              Q2_up = quantile(Q2, 0.975),
              Q3_low = quantile(Q3, 0.025),
              Q3_med = quantile(Q3, 0.5),
              Q3_up = quantile(Q3, 0.975)
    ) %>%
    mutate(mod = Model)
  
  med_data_plot <- adastra_dat_analysis %>%
    group_by(Timepoint_ID) %>%
    summarise(
      q1 = quantile(daily_VL, 0.25),
      q2 = quantile(daily_VL, 0.5),
      q3 = quantile(daily_VL, 0.75)
    ) %>%
    mutate(mod = Model)
  
  
  list(adastra_dat_analysis, preds3, med_data_plot)
  
}

list_data_expo <- posterior_check_plot(model_settings, "Exponential", Dmax = 5.5, Type = "All")
list_data_biexpo <- posterior_check_plot(model_settings, "Bi_exponential", Dmax = 5.5, Type = "All")
list_data_growth_decay <- posterior_check_plot(model_settings, "Up_and_down", Dmax = 5.5, Type = "All")

adastra_dat_analysis <- rbind(list_data_expo[[1]], list_data_biexpo[[1]],list_data_growth_decay[[1]])
preds3 <- rbind(list_data_expo[[2]], list_data_biexpo[[2]],list_data_growth_decay[[2]])
med_data_plot <- rbind(list_data_expo[[3]], list_data_biexpo[[3]],list_data_growth_decay[[3]])

adastra_dat_analysis$mod <- as.factor(adastra_dat_analysis$mod)
adastra_dat_analysis$mod <- factor(adastra_dat_analysis$mod, levels = levels(adastra_dat_analysis$mod)[c(2,1,3)])
levels(adastra_dat_analysis$mod) <- c("Exponential decay model",
                                      "Bi-exponential decay model",
                                      "Growth-and-decay model")

preds3$mod <- as.factor(preds3$mod)
preds3$mod <- factor(preds3$mod, levels = levels(preds3$mod)[c(2,1,3)])
levels(preds3$mod) <- c("Exponential decay model",
                        "Bi-exponential decay model",
                        "Growth-and-decay model")

med_data_plot$mod <- as.factor(med_data_plot$mod)
med_data_plot$mod <- factor(med_data_plot$mod, levels = levels(med_data_plot$mod)[c(2,1,3)])
levels(med_data_plot$mod) <- c("Exponential decay model",
                               "Bi-exponential decay model",
                               "Growth-and-decay model")


G1 <- ggplot(mapping = aes(x = Timepoint_ID)) +
  geom_hline(yintercept = 0, col = "red", linetype = "21", linewidth = 0.5) +
  geom_jitter(data = adastra_dat_analysis, col = "gray40",
              aes(y = daily_VL, shape = censor), alpha = 0.3, size = 2,
              position = position_jitter(width = 0.2, seed = 1)) +
  scale_shape_manual(values = c(6, 1), guide = NULL) +
  #Median
  geom_ribbon(data = preds3, aes(ymin = Q2_low  , ymax = Q2_up, fill = mod), alpha = 0.35) +
  geom_line(data = preds3,aes(y = Q2_med ,col = mod), linewidth = 1.3) +
  #Q1
  geom_ribbon(data = preds3, aes(ymin = Q1_low, ymax = Q1_up, fill = mod), alpha = 0.35) +
  geom_line(data = preds3, aes(y = Q1_med, col = mod), linetype = '31',linewidth = 0.75) +
  #Q3
  geom_ribbon(data = preds3, aes(ymin = Q3_low, ymax = Q3_up, fill = mod), alpha = 0.35) +
  geom_line(data = preds3, aes(y = Q3_med, col = mod), linetype = '31',linewidth = 0.75) +
  #Summary observed
  geom_point(data = med_data_plot, aes(y = q2), shape = 17, size = 3.5,
             col = 'gray20') +
  geom_errorbar(data = med_data_plot, aes(ymin = q1, ymax = q3),
                width = 0.2, linewidth = 0.5) +
  theme_bw(base_size = 15) +
  coord_cartesian(ylim = c(0, 8)) +
  scale_x_continuous(breaks = seq(0, 5, by = 1)) +
  ylab("Viral densities (log10 genomes/mL)") +
  xlab("Time since randomisation (days)") +
  ggtitle("A) Predicted clearance kinetics") +
  theme(axis.title = element_text(face = 'bold')) +
  facet_grid(.~mod) +
  scale_color_manual(values = c("#0A5EB0", "#F72C5B", "#F26B0F"), guide = NULL) +
  scale_fill_manual(values = c("#0A5EB0", "#F72C5B", "#F26B0F"), guide = NULL) 
G1



adastra_dat_analysis <- adastra_dat_analysis %>%
  group_by(ID, Timepoint_ID, mod) %>%
  mutate(rediduals = log10_viral_load - preds) %>%
  ungroup()

G2 <- ggplot(adastra_dat_analysis %>% filter(censor == "none"), aes(x = as.factor(Timepoint_ID), y = rediduals)) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 2, aes(shape = censor, fill = censor, col=censor)
  ) +
  geom_boxplot(outlier.shape = NA, fill = "white", alpha = 0.6) +
  theme_bw(base_size = 15) +
  facet_grid(.~ mod) +
  geom_hline(yintercept = 0, col = "red", linetype = "21", linewidth = 0.5) +
  xlab("Time since randomisation (days)") +
  ylab("Residuals") +
  theme(
    axis.title  = element_text(face = "bold")
  ) +
  scale_shape_manual(values = c(21), name = "Censor", guide = NULL) +
  scale_fill_manual(values = rev(c("#38419D")), name = "Censor", guide = NULL) +
  scale_color_manual(values = rev(c("#38419D")), name = "Censor", guide = NULL) +
  ggtitle("B) Residuals")

G2



png(here("Plots", "R1", "prediction_residuals.png"), width = 10, height = 10, units = 'in', res = 350)
ggarrange(G1, G2, ncol = 1)
dev.off()  
