model_settings <- read.csv("Outputs/model_settings.csv")
Dmax <- 5.5
Type <- "All"

ind <- which(model_settings$Dmax == Dmax & model_settings$Type == Type & grepl("Up_and_down", model_settings$mod))
load(paste0("Rout/model_settings_job_", ind, ".RData"))
load(paste0("Rout/model_fits_job_", ind, ".RData"))

thetas <- rstan::extract(out)
tmax <- thetas$theta_rand_id[,,3] %>% apply(2, function (x) x + thetas$tmax_pop) %>%
  apply(2, median)

sum(tmax < 0)/(length(tmax)) 
sum(tmax > 1)/(length(tmax)) 

####################################################################################################
model_settings <- read.csv("Outputs/model_settings.csv")
Dmax <- 5.5
Type <- "All"

ind <- which(model_settings$Dmax == Dmax & model_settings$Type == Type & grepl("Bi_exponential", model_settings$mod))
load(paste0("Rout/model_settings_job_", ind, ".RData"))
load(paste0("Rout/model_fits_job_", ind, ".RData"))

thetas <- rstan::extract(out)

intercept_1 <- thetas$theta_rand_id[,,3] %>% apply(2, function (x) x + thetas$intercept[,1]) %>%
  apply(2, median)

sum(intercept_1 < 2)/(length(intercept_1)) 

hist(intercept_1)


sum(tmax < 0)/(length(tmax)) 
sum(tmax > 1)/(length(tmax)) 