library(tidyverse)
library(cmstanr)
library(loo)
library(bayesplot)

setwd("path/to/modeling-of-control")

#load model 1
independent_fit <- read_rds("model_fitting/2D_model_replication/fit.rds")
independent_loo <- loo::loo(independent_fit$draws("log_lik"), cores=4, save_psis=TRUE)
print(independent_loo)
plot(independent_loo)

#load model 2
tradeoff_fit <- read_rds("model_fitting/1D_model_replication/fit.rds")
tradeoff_loo <- loo::loo(tradeoff_fit$draws("log_lik"), cores=4, save_psis=TRUE)
print(tradeoff_loo)
plot(tradeoff_loo)

#compare model fits
loo::loo_compare(independent_loo, tradeoff_loo)


