library(tidyverse)
library(tidybayes)
library(bayestestR)
library(bayesplot)
library(RWiener)
library(purrr)
library(cmdstanr)

# load model fit
setwd("path/to/modeling-of-control")
path <- "model_fitting/NoTradeoff_model/"
fit <- read_rds(paste0(path, "fit.rds"))

df <- read_csv("data/dataset1.csv", col_types = cols()) %>% 
  select(trialCount, block, blockType, acc, RT, stimCongruency, switchType, subject) %>% 
  mutate(RT = RT/1000,
         incongruency = ifelse(stimCongruency == "i", 1, 0),
         taskSequence = ifelse(switchType == "s", 1, 0),
         switchProp = ifelse(blockType == "B" | blockType == "D", 0.75, 0.25),
         incProp = ifelse(blockType == "A" | blockType == "B", 0.75, 0.25),
         subject = match(subject, unique(subject))) %>% 
  filter(RT < 1.5 & RT > 0.3)

conditions <- df %>% 
  distinct(incongruency, taskSequence, switchProp, incProp)

conditions <- model.matrix(~ 0 + incongruency + taskSequence + switchProp + incProp +
                             incongruency:taskSequence + incongruency:incProp + taskSequence:switchProp, conditions) 

desired_order <- c(
  "incongruency",
  "taskSequence",
  "switchProp",
  "incProp",
  "incongruency:taskSequence",
  "incongruency:incProp",
  "taskSequence:switchProp"
)
conditions <- conditions[, desired_order]

stan_data <- list(
  N = nrow(df),
  RTs = df$RT,
  acc = df$acc,
  incongruencies = df$incongruency,
  switches = df$taskSequence,
  K = 7,
  blockN = df$block,
  nSubjects = length(unique(df$subject)),
  subjects = df$subject,
  priorOnly = 0,
  X_pred = conditions,
  N_pred = nrow(conditions)
)

# fit model to get draws 
model <- cmdstan_model(stan_file = paste0(path, "model_pred.stan"))
fit.pred <- model$generate_quantities(fit, data=stan_data)
draws <- spread_draws(fit.pred, final_drift[N], final_bs[N], final_ndt[N], final_bias[N])

# expected reaction time
posterior_epred_wiener <- function(bs, ndt, bias, mu) {
  # obtained from https://doi.org/10.1016/j.jmp.2009.01.006
  # mu is the drift rate
  ndt - (bias * bs) / mu + bs / mu *
    (exp(-2 * mu * (bias * bs)) - 1) / (exp(-2 * mu * bs) - 1)
}

posterior_pred <- draws %>%
  mutate(RT = pmap_dbl(list(final_bs, final_ndt,
                            final_bias, final_drift), posterior_epred_wiener))

#add in details about condition for each row
task_conditions <- as_tibble(conditions) %>% mutate(N=row_number())

# final posterior pred data frame, with cols for switchProp/incProp
posterior_pred <- posterior_pred %>% 
  left_join(task_conditions, by="N") %>% 
  mutate(congruency = ifelse(incongruency == 1, "i", "c"),
         taskSequence = ifelse(taskSequence == 1, "s", "r"),
         switchProp = ifelse(switchProp == 0.75, "75%", "25%"),
         incProp = ifelse(incProp == 0.75, "75%", "25%")) %>% 
  select(N, .chain, .iteration, .draw, congruency, taskSequence, switchProp, incProp, RT)

write_csv(posterior_pred, paste0(path, "posterior_pred.csv"))
