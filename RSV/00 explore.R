day <- 1:6
mean_placebo <- c(7.5, 6.8, 6.3, 6, 5.6, 5.2)
mean_drug <- c(7.4, 6.6, 5.8, 5.2, 4.9, 4.6)

sd <- 1.2

n <- 100
n_rep <- 2


mean <- mean_placebo
group <- "placebo"


sim_pop <- function(day, mean, sd, n, n_rep, group){
  output <- NULL
  for(i in day){
    temp <- data.frame(
      "day" = i,
      "log10_vl" = rnorm(n, mean = mean[i], sd),
      "arm" = group )      
    output <- rbind(output, temp)
  }
  output
}


sims <- rbind(
sim_pop(day, mean_placebo, sd, n, n_rep, "Placebo"),
sim_pop(day, mean_drug, sd, n, n_rep, "Drug")
)

sims_median <- sims %>% group_by(arm, day) %>% summarise(median = median(log10_vl),
                                                         mean = mean(log10_vl))


ggplot(sims) +
  geom_jitter(aes(x = day, y = log10_vl, col = arm), width = 0.2, alpha = 0.25, size = 2) +
  theme_bw(base_size = 13) +
  ylab("Viral density (log10 genomes/mL)") +
  xlab("Time point (days)") +
  scale_x_continuous(breaks = 0:8) +
  scale_y_continuous(breaks = seq(0,12,2)) +
  scale_color_manual(values = c("#0D92F4", "#AF1740"), name = "") +
  geom_line(data = sims_median, aes(x = day , y = mean, col = arm), linetype = "dashed", linewidth = 0.6) +
  geom_point(data = sims_median, aes(x = day , y = mean, col = arm), shape = 17, size = 4) +
  geom_hline(yintercept = 2, linetype = "dashed", col = "red", linewidth = 0.75)
