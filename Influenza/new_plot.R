Model <- "Exponential"
color <- "#0A5EB0"


posterior_check_plot <- function(Model, color){
  ind_exp <- which(model_settings$Dmax == Dmax & model_settings$Type == Type &
                     grepl(Model, model_settings$mod))
  
  load(paste0("Rout/model_settings_job_", ind_exp, ".RData"))
  load(paste0("Rout/model_fits_job_", ind_exp, ".RData"))
  
  data_plot[[1]] <- adastra_dat_analysis
  
  preds <-  rstan::extract(out, pars = 'preds')
  data_plot[[1]]$preds <- apply(preds$preds, 2, median) 
  
  preds2 <- cbind("ID" = data_plot[[1]]$ID,
                  "Timepoint_ID" = data_plot[[1]]$Timepoint_ID, 
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
    )
  
  med_data_plot <- data_plot[[1]] %>%
    group_by(Timepoint_ID) %>%
    summarise(
      q1 = quantile(daily_VL, 0.25),
      q2 = quantile(daily_VL, 0.5),
      q3 = quantile(daily_VL, 0.75)
    )


  G <- ggplot(mapping = aes(x = Timepoint_ID)) +
    geom_hline(yintercept = 0, col = 'red', linetype = 'dashed',
               linewidth = 0.75) +
    geom_jitter(data = data_plot[[1]], col = "gray30",
                aes(y = daily_VL, shape = censor), alpha = 0.3, size = 2,
                position = position_jitter(width = 0.2, seed = 1)) +
    scale_shape_manual(values = c(6, 1), guide = NULL) +
    geom_point(data = med_data_plot, aes(y = q2), shape = 17, size = 3,
               col = 'gray20') +
    geom_errorbar(data = med_data_plot, aes(ymin = q1, ymax = q3),
                  width = 0.2, linewidth = 0.75) +
    geom_line(data = med_data_plot,
              aes(y = q2), col = "gray30", linetype = 'dashed',
              linewidth = 0.75) +
    geom_ribbon(data = preds3, fill = color, col = color, linetype = 'dashed',
      aes(ymin = Q1_med  , ymax = Q3_med), alpha = 0.2) +
    geom_line(data = preds3, col = color,
              aes(y = Q2_med), alpha = 0.75, linewidth = 2) +
    theme_bw(base_size = 13) +
    coord_cartesian(ylim = c(0, 8)) +
    scale_x_continuous(breaks = seq(0, 5, by = 1)) +
    xlab("") +
    ylab("") +
    theme(title = element_text(size = 12))
    
    G
  
}

G1 <- posterior_check_plot("Exponential", "#0A5EB0") + ggtitle("Exponential decay")
G2 <- posterior_check_plot("Bi_exponential", "#F72C5B") + ggtitle("Bi-exponential decay")
G3 <- posterior_check_plot("Up_and_down", "#F26B0F") + ggtitle("Growth-and-decay")

library(ggpubr)

annotate_figure(
  ggarrange(G1, G2, G3, common.legend = T, legend = "bottom", ncol = 3),
  left = text_grob("Viral densities \n(log10 genomes/mL)", rot = 90, vjust = 1, size = 13,
                   face = 'bold'),
  bottom = text_grob("Time since randomisation (days)", size = 13, face = 'bold'),
  top = text_grob("A) Predicted clearance kinetics", face = "bold", size = 16,
                  hjust = 1)
  
)


