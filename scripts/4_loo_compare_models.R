library(tidyverse)
library(cmdstanr)
library(loo)
library(bayesplot)

setwd("path/to/modeling-of-control")

#load model 1
independent_fit <- read_rds("model_fitting/2D_model/fit.rds")
independent_loo <- loo::loo(independent_fit$draws("log_lik"), cores=4, save_psis=TRUE)
print(independent_loo)
plot(independent_loo)

#load model 2
tradeoff_fit <- read_rds("model_fitting/1D_model/fit.rds")
tradeoff_loo <- loo::loo(tradeoff_fit$draws("log_lik"), cores=4, save_psis=TRUE)
print(tradeoff_loo)
plot(tradeoff_loo)

#compare model fits
loo::loo_compare(independent_loo, tradeoff_loo)


#-------------------
# Analysis with 3rd forced-no-tradeoff model
model3_fit <- read_rds("model_fitting/NoTradeoff_model/fit.rds")
model3_loo <- loo::loo(model3_fit$draws("log_lik"), cores=4, save_psis=TRUE)
print(model3_loo)
plot(model3_loo)

#compare model fits
loo::loo_compare(independent_loo, model3_loo)
loo::loo_compare(independent_loo, tradeoff_loo, model3_loo)
