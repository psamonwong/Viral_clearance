library(ggpubr)

bi_exp <- function(a0, b0, alpha, beta, tmax){
  time = seq(0,tmax,0.1)
  
  res =  data.frame(time = time,
                    y = log(exp(a0 - alpha*time) + exp(b0 - beta*time)),
                    y1 = a0 - alpha * time,
                    y2 = b0 - beta * time
  )
  
  G <- ggplot(res, aes(x = time)) +
    geom_line(aes(y = y), linewidth = 1) +
    geom_line(aes(y = y1), linewidth = 0.5, col = "red", linetype = "dashed") +
    geom_line(aes(y = y2), linewidth = 0.5, col = "blue", linetype = "dashed") +
    ylim(0, 8) +
    theme_bw(base_size = 13) +
    scale_x_continuous(breaks = 0:tmax, limits = c(0,tmax))
  
  return(G)
  
}


bi_exp(a0 = 6, b0 = 2, alpha = 1, beta = 0.1, tmax = 14)


ggarrange(bi_exp(a0 = 6, b0 = 2, alpha = 1, beta = 0.1, tmax = 14),
          bi_exp(a0 = 6, b0 = 2, alpha = 1*2, beta = 0.1, tmax = 14),
          bi_exp(a0 = 6, b0 = 2, alpha = 1, beta = 0.1*2, tmax = 14),
          bi_exp(a0 = 6, b0 = 2, alpha = 1*2, beta = 0.1*2, tmax = 14),
          ncol = 2, nrow = 2, labels = "AUTO", align = "hv")
