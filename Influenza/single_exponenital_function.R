library(ggpubr)

single_exp <- function(a0, alpha, tmax){
  time = seq(0,tmax,0.1)
  
  res =  data.frame(time = time,
                    y = a0 - alpha * time
  )
  
  G <- ggplot(res, aes(x = time)) +
    geom_line(aes(y = y), linewidth = 1) +
  #  geom_line(aes(y = y1), linewidth = 0.5, col = "red", linetype = "dashed") +
  # geom_line(aes(y = y2), linewidth = 0.5, col = "blue", linetype = "dashed") +
    ylim(0, 8) +
    theme_bw(base_size = 13) +
    scale_x_continuous(breaks = seq(0,tmax, 2), limits = c(0,tmax)) +
    xlab("Time post-randomisation (days)") +
    ylab("Log viral densities/mL")
  
  return(G)
  
}


png("Plots/single_exponential_example.png", width = 4, height = 4, units = "in", res = 350)
single_exp(a0 = 6, alpha = 1, tmax = 14)
dev.off()
