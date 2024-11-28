ff = list.files('Rout/', pattern = 'A_and_B_5d',)
ff = ff[grep(pattern = 'model_fits_',x = ff, ignore.case = T)]
if(!length(ff)==nrow(model_settings)) stop('not all outputs are ready for all model settings')
ff = paste0('Rout/',ff)


model_list = list()
mod_select <- c(1,2,3)

stan_inputs_i <- stan_inputs
platcov_dat_analysis <- data.frame("ID" = stan_inputs_i$analysis_data_stan$id,
                                   "log10_viral_load" = stan_inputs_i$analysis_data_stan$log_10_vl,
                                   "Time" = stan_inputs_i$analysis_data_stan$obs_day)

ID_map <- stan_inputs_i$ID_map

for(j in 1:length(mod_select)){
  load(ff[j])
  model_list[[j]] = out
}

preds <- lapply(model_list, rstan::extract, "preds")
preds_list <- list()
mods <- c("Linear", "Bi-exponential", "Up and down")

for(j in 1:length(mod_select)){
  preds_list[[j]] <- preds[[j]]$preds %>% apply(2, quantile, c(0.025, 0.5, 0.975)) %>% t()
  preds_list[[j]] <- cbind(platcov_dat_analysis,  preds_list[[j]])
  preds_list[[j]]$model <- mods[j]
  colnames(preds_list[[j]]) <- c(colnames(platcov_dat_analysis), c("low", "med", "up", "model"))
  preds_list[[j]]$censor <- "left"
  preds_list[[j]]$censor[1:stan_inputs_i$analysis_data_stan$N_obs] <- "None"
  preds_list[[j]] <- merge(preds_list[[j]], ID_map, by.x = "ID", by.y = "ID_stan")
  preds_list[[j]] <- merge(preds_list[[j]], as.data.frame(Baseline_data[,c("ID", "Trt")]), by.x = "ID_key", by.y = "ID")

}


preds_dat  <- do.call("rbind", preds_list)
preds_dat$model <- as.factor(preds_dat$model)
preds_dat$censor <- as.factor(preds_dat$censor)

ID_map = merge(stan_inputs_i$ID_map, Baseline_data, by.x = 'ID_key',by.y = 'ID')
ID_map <- ID_map[order(ID_map$fluType, ID_map$Trt, ID_map$ID_key),]


ind_plot_list <- list()
resid_dat <- NULL
for(i in 1:nrow(ID_map)){
  plot_data <- preds_dat %>% filter(ID_key == ID_map$ID_key[i])
  plot_data$resid <- plot_data$log10_viral_load - plot_data$med
  plot_data$Timepoint_ID <- round(plot_data$Time)
  resid_dat <- rbind(resid_dat, plot_data)
  
  lab <- paste0(plot_data$ID_key[1], "\nFlu ", ID_map$fluType[i])
  
  ind_plot_list[[i]] <- ggplot() +
    geom_point(data = plot_data, aes(x = Time, y = log10_viral_load, shape = censor),
               size = 2.5, alpha = 0.7) +
    geom_ribbon(data = plot_data, aes(x = Time, ymin = low, ymax = up, fill = model), alpha = 0.2) +
    geom_line(data = plot_data, aes(x = Time, y = med, col = model), linewidth = 0.75) +
    theme_bw() +
    scale_y_continuous(labels=label_math(), breaks = seq(0,10,2)) +
    coord_cartesian(ylim = c(0,9), xlim = c(0,5))+
    scale_x_continuous(breaks = 0:14) +
    ylab("") +
    xlab("") +
    theme(
      axis.title  = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 8),
      legend.position = "right",
      plot.margin = unit(c(0.1,0.1,0.1,0.1), 'lines')) +
    scale_color_manual(values = c("#1640D6", "#BE3144", "black"), name = "Model") +
    scale_fill_manual(values = c("#1640D6", "#BE3144",  "black"), name = "Model") +
    scale_shape_manual(values = c(17, 16), guide = "none", drop=FALSE) +
    ggtitle(lab)
  
}

ind_plot_all <- ggarrange(plotlist =  ind_plot_list, nrow = 4, ncol = 4, common.legend = T, legend = "right")

for(i in 1:length(ind_plot_all)){
  png(paste0("Plots/Individual_plots/Ind_plot_", i, ".png"), width = 12, height = 8, units = "in", res = 350)
  print(annotate_figure( ind_plot_all[[i]], bottom = textGrob("Time since randomisation (days)", vjust = 0.5, gp = gpar(cex = 1.2, fontface="bold")),
                         left = textGrob("Viral genomes/mL", rot = 90, gp = gpar(cex = 1.2, fontface="bold"))))
  dev.off()
}


resid_dat$model <- factor(resid_dat$model, levels = c("Linear", "Bi-exponential", "Up and down"))

resid_plot <- ggplot(resid_dat[resid_dat$Timepoint_ID <= 5,], aes(x = as.factor(Timepoint_ID), y = resid)) +
  geom_jitter(width = 0.2, alpha = 0.45, size = 2, aes(shape = censor, color = censor)) +
  geom_boxplot(outlier.shape = NA, fill = "white", alpha = 0.6) +
  theme_bw(base_size = 14) +
  facet_grid(.~ model) +
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
  xlab("Time since randomisation (days)") +
  ylab("Residuals") +
  theme(
   # axis.title  = element_text(face = "bold", size = 12),
  #  strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom"
  ) +
  scale_shape_manual(values = c(17, 16), name = "Censor") +
  scale_color_manual(values = rev(c("#38419D", "#BF3131")), name = "Censor")
resid_plot 

png("Plots/residual_plot.png" , width = 8, height = 6, units = "in", res = 350)
resid_plot
dev.off()


