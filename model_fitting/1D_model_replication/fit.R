library(cmdstanr)
library(tidyverse)
library(tidybayes)
library(bayesplot)

setwd("path/to/modeling-of-control")

# rstan options
options(mc.cores = parallel::detectCores())

#load data
df <- read_csv("data/dataset2_replication.csv", col_types = cols()) %>% 
  mutate(RT = RT/1000,
         congruency = ifelse(congruency == "i", 1, 0),
         taskSequence = ifelse(taskSequence == "s", 1, 0),
         subject = match(subject, unique(subject))) %>% 
  filter(RT < 1.5 & RT > 0.3)

stan_data <- list(
  N = nrow(df),
  RTs = df$RT,
  acc = df$acc,
  incongruencies = df$congruency,
  switches = df$taskSequence,
  K = 7,
  blockN = df$block,
  nSubjects = length(unique(df$subject)),
  subjects = df$subject,
  priorOnly = 0
)

path <- "model_fitting/1D_model_replication/"
model <- cmdstan_model(stan_file = paste0(path, "model.stan"))

# initialize best on optimizing best parameter values
optim_init_fun <- function(){
  list(
    Intercept_d = 2,
    Intercept_bs = 0.7,
    Intercept_ndt = -4,
    b_d = rep(0, 7),
    b_bs = rep(0, 7),
    b_ndt = rep(0, 7),
    sd_d = 0.001,
    sd_bs = 0.001,
    sd_ndt = 0.001,
    z_d = rep(0.001, length(unique(df$subject))),
    z_bs = rep(0.001, length(unique(df$subject))),
    z_ndt = rep(0.001, length(unique(df$subject)))
  )
}
optim_fit <- model$optimize(data=stan_data, init=optim_init_fun, iter=4000)
optim_params_df <- optim_fit$summary() 
optim_params <- list(
  Intercept_d = optim_params_df[optim_params_df$variable == "Intercept_d", ]$estimate,
  Intercept_bs = optim_params_df[optim_params_df$variable == "Intercept_bs", ]$estimate,
  Intercept_ndt = optim_params_df[optim_params_df$variable == "Intercept_ndt", ]$estimate,
  b_d = optim_params_df[optim_params_df$variable %in% rep(paste0("b_d[", 1:7, "]"), each=1), ]$estimate,
  b_bs = optim_params_df[optim_params_df$variable %in% rep(paste0("b_bs[", 1:7, "]"), each=1), ]$estimate,
  b_ndt = optim_params_df[optim_params_df$variable %in% rep(paste0("b_ndt[", 1:7, "]"), each=1), ]$estimate,
  sd_d = optim_params_df[optim_params_df$variable == "sd_d", ]$estimate,
  sd_bs = optim_params_df[optim_params_df$variable == "sd_bs", ]$estimate,
  sd_ndt = optim_params_df[optim_params_df$variable == "sd_ndt", ]$estimate,
  z_d = optim_params_df[optim_params_df$variable %in% rep(paste0("z_d[", 1:length(unique(df$subject)), "]"), each=1), ]$estimate,
  z_bs = optim_params_df[optim_params_df$variable %in% rep(paste0("z_bs[", 1:length(unique(df$subject)), "]"), each=1), ]$estimate,
  z_ndt = optim_params_df[optim_params_df$variable %in% rep(paste0("z_ndt[", 1:length(unique(df$subject)), "]"), each=1), ]$estimate
)

save(optim_params, file=paste0(path, "optim_params_for_inits.rda"))
optim_fit$save_object(file = paste0(path, "optim_fit.rds"))

# OR

load(paste0(path, "optim_params_for_inits.rda"))

initfun <- function() {optim_params}

fit <- model$sample(data=stan_data, chains=4, iter_warmup=500, init=initfun, iter_sampling=1000, max_treedepth=12, refresh=10)
fit$save_object(file = paste0(path, "fit.rds"))