library(cmdstanr)
library(tidyverse)
library(tidybayes)
library(bayesplot)

#load model
setwd("path/to/modeling-of-control")
path <- "model_fitting/1D_model/"
fit <- read_rds(paste0(path, "fit.rds"))

# get model summaries
dignostic_summary <- fit$diagnostic_summary()
model_summary <- fit$summary()

# print 5 largest RHats
model_summary %>% 
  arrange(desc(rhat)) %>%
  slice_head(n = 5)

#print 5 smallest bulkess
model_summary %>% 
  arrange(ess_bulk) %>%
  slice_head(n = 5)

#traceplots
mcmc_trace(fit$draws(paste0("Intercept_", c("d", "bs", "ndt"))), facet_args = list(nrow = 3))
mcmc_trace(fit$draws(rep(paste0("b_d[", 1:7, "]"))))
mcmc_trace(fit$draws(rep(paste0("b_bs[", 1:7, "]"), each=1)))
mcmc_trace(fit$draws(rep(paste0("b_ndt[", 1:7, "]"), each=1)))

#hists
mcmc_hist(fit$draws(paste0("Intercept_", c("d", "bs", "ndt"))))
mcmc_hist(fit$draws(rep(paste0("b_d[", 1:15, "]"), each=1)))
mcmc_hist(fit$draws(rep(paste0("b_bs[", 1:15, "]"), each=1)))
mcmc_hist(fit$draws(rep(paste0("b_ndt[", 1:15, "]"), each=1)))
mcmc_hist(fit$draws(rep(paste0("z_d[", 1:60, "]"))))
mcmc_hist(fit$draws(rep(paste0("z_ndt[", 1:60, "]"))))
mcmc_hist(fit$draws(rep(paste0("z_bs[", 1:60, "]"))))

# pairs
mcmc_pairs(fit$draws(), pars = c("b_d[1]", "b_bs[1]", 'b_ndt[1]'))
mcmc_pairs(fit$draws(c("b_d[2]", "b_bs[2]", 'b_ndt[2]')))

# make plot of all main drift 
drift_vars <- c("Intercept_d", rep(paste0("b_d[", 1:7, "]")))
bs_vars <- c("Intercept_bs", rep(paste0("b_bs[", 1:7, "]")))
ndt_vars <- c("Intercept_ndt", rep(paste0("b_ndt[", 1:7, "]")))

all_vars <- c(drift_vars, bs_vars, ndt_vars)
all_vars_matrix <- matrix(all_vars, nrow = 3, byrow = TRUE)
all_vars_matrix_inverted <- aperm(all_vars_matrix, c(2, 1))
inverted_vars_list <- unlist(apply(all_vars_matrix_inverted, 1, list))

hist_plot <- mcmc_trace(fit$draws(inverted_vars_list), facet_args=list(ncol=3)) +
  ggtitle("1D Model Parameter Traceplots") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())
hist_plot

ggsave("plots/1d_traceplot.jpg", hist_plot, dpi = 300, width = 6.5, height = 8, units = "in")

