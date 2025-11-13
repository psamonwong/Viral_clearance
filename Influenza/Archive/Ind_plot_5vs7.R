model_settings <- read.csv("Outputs/model_settings.csv")
Dmax <- c(5.5, 7.5)
Type <- "All"

ind <- which(model_settings$Dmax %in% Dmax & model_settings$Type == Type & grepl("Exponential_decay", model_settings$mod))

data_plot <- list()

for(i in 1:length(ind)){
  load(paste0("Rout/model_settings_job_", ind[i], ".RData"))
  load(paste0("Rout/model_fits_job_", ind[i], ".RData"))
  
  data_plot[[i]] <- adastra_dat_analysis
  
  preds <-  rstan::extract(out, pars = 'preds')
  data_plot[[i]]$preds <- apply(preds$preds, 2, median)
  data_plot[[i]]$Dmax <- Dmax[i]
}

data_Plot <- do.call("rbind", data_plot)

#Individual plots
data_Plot$censor <- as.factor(data_Plot$censor)
data_Plot$Dmax <- as.factor(floor(data_Plot$Dmax))
levels(data_Plot$Dmax) <- paste0(levels(data_Plot$Dmax), " days")


IDs <-data_Plot %>% 
  arrange(fluType, ID, Time) %>%
  pull(ID) %>%
  unique()

ind_plot_list <- list()

for(i in 1:length(IDs)){
  data_Plot_i <- data_Plot %>% filter(ID == IDs[i])
  data_Plot_median_i <- data_Plot_i %>% group_by(Timepoint_ID) %>%
    mutate(med_vl = median(log10_viral_load)) %>%
    ungroup() %>%
    distinct(Timepoint_ID, .keep_all = T) %>%
    select(Time, med_vl)
  
  lab <- paste0(data_Plot_i$ID[1], "\nInfluenza ", data_Plot_i$fluType[1])
  
  ind_plot_list[[i]] <- 
    ggplot(data_Plot_i, aes(x = Time)) +
    geom_point(aes(y = log10_viral_load, shape = censor), size = 1.75) +
    geom_line(aes(y = preds, col = Dmax), linewidth = 1, alpha = 0.75) +
    #geom_point(data = data_Plot_median_i, mapping = aes(x = Time, y = med_vl),
    #           size=1.75, shape = 17, alpha = 0.75) +
    #geom_line(data = data_Plot_median_i, mapping = aes(x = Time, y = med_vl),
    #          linetype = "dashed", , alpha = 0.75) +
    scale_color_manual(values = c("#69247C", "#DA498D","#F29F58"), name = "Follow-up") +
    theme_bw(base_size = 10)   +  
    scale_y_continuous(labels=label_math(), breaks = seq(0,8,2)) +
    coord_cartesian(ylim = c(0,8), xlim = c(0,7)) +
    scale_shape_manual(values = c(25, 21), guide = "none", drop=FALSE) +
    geom_hline(yintercept = 0, col = 'red', linetype = "dashed") +
    ylab("") +
    xlab("")  +
    ggtitle(lab) +
    theme(plot.title = element_text(face = "bold",hjust = 0.5, size = 8))
  
}

ind_plot_all <- ggarrange(plotlist =  ind_plot_list, nrow = 7, ncol = 4, common.legend = T, legend = "right")

plot_now = T

if(plot_now){
  for(i in 1:length(ind_plot_all)){
    png(paste0("Plots/Individual_plots_5vs7/Ind_plot_", i, ".png"), width = 12, height = 12, units = "in", res = 350)
    print(annotate_figure( ind_plot_all[[i]], bottom = textGrob("Time since randomisation (days)", vjust = 0.5, gp = gpar(cex = 1.2, fontface="bold")),
                           left = textGrob("Viral densities (genomes/mL)", rot = 90, gp = gpar(cex = 1.2, fontface="bold"))))
    dev.off()
  }  
}




