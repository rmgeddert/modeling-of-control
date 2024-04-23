library(cmdstanr)
library(tidyverse)
library(tidybayes)
library(bayesplot)
library(bayestestR)
library(patchwork)

setwd("path/to/modeling-of-control")

posterior_plot_func <- function(posterior_df, y_par, color_labels){
  pd <- position_dodge(0.2)
  ggplot(posterior_df, aes(x=taskSequence, y=!!rlang::sym(y_par), color=congruency, group=congruency)) +
    geom_pointinterval(aes(ymin=.lower, ymax=.upper), position=pd) +
    geom_line(position=pd, show.legend = FALSE) +
    facet_grid(switchProp ~ incProp,
               labeller = labeller(incProp = c("25%" = "25% Incongruent",
                                               "75%" = "75% Incongruent"),
                                   switchProp = c("25%" = "25% Switch",
                                                  "75%" = "75% Switch"))) +
    scale_color_manual(name="Congruency",
                       values = color_labels,
                       labels = c("c"="Congruent", "i"="Incongruent")) +
    scale_x_discrete(labels=c("r" = "Repeat", "s" = "Switch")) +
    xlab("Task Sequence") +
    guides(color = guide_legend(override.aes = list(point_size=4,
                                                    linetype="blank"))) +
    theme(axis.title.y = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5)) 
}

effects_plot_func <- function(effect_df){
  ggplot(effect_df, aes(y = effect, x = .value, xmin=.lower, xmax = .upper, color=as.factor(sig))) +
    geom_pointinterval(position=position_dodge(width=0.7)) +
    geom_vline(xintercept = 0, linetype="dashed") +
    scale_y_discrete(
      labels = c(
        "TaskSequence" = "Task Sequence",
        "Congruency" = "Congruency",
        "TaskSequence:SwitchProp" = "Task Sequence x\n Switch Prop",
        "TaskSequence:IncProp" = "Task Sequence x\n Inc Prop",
        "Congruency:SwitchProp" = "Congruency x\n Switch Prop",
        "Congruency:IncProp" = "Congruency x\n Inc Prop"
      )
    ) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none") +
    guides(color = guide_legend(override.aes = list(size=10))) +
    scale_color_manual(values = c("#56B1F7", "#132B43"),
                         breaks=c(1, 0))
}

rename_list = c(
  "Intercept" = "(Intercept)",
  "Congruency" = "congruencyi",
  "TaskSequence" = "taskSequences",
  "SwitchProp" = "switchProp75%",
  "IncProp" = "incProp75%",
  "Congruency:TaskSequence" = "congruencyi:taskSequences",
  "Congruency:SwitchProp" = "congruencyi:switchProp75%",
  "TaskSequence:SwitchProp" = "taskSequences:switchProp75%",
  "Congruency:IncProp" = "congruencyi:incProp75%",
  "TaskSequence:IncProp" = "taskSequences:incProp75%",
  "SwitchProp:IncProp" = "switchProp75%:incProp75%",
  "Congruency:TaskSequence:SwitchProp" = "congruencyi:taskSequences:switchProp75%",
  "Congruency:TaskSequence:IncProp" = "congruencyi:taskSequences:incProp75%",
  "Congruency:SwitchProp:IncProp" = "congruencyi:switchProp75%:incProp75%",
  "TaskSequence:SwitchProp:IncProp" = "taskSequences:switchProp75%:incProp75%",
  "Congruency:TaskSequence:SwitchProp:IncProp" = "congruencyi:taskSequences:switchProp75%:incProp75%"
)

include_in_plot_list = c(
  "TaskSequence",
  "Congruency",
  "TaskSequence:SwitchProp",
  "TaskSequence:IncProp",
  "Congruency:SwitchProp",
  "Congruency:IncProp"
)

#--------------------------------------------------
# Load things
path = "model_fitting/2D_model_replication/"

# load model DDM posteriors
DDM_posteriors <- read_csv(paste0(path, "DDM_posteriors.csv"))

#--------------------------------------------------
# Get Median QIs for each DDM par
#DRIFT
drift_qi <- DDM_posteriors %>% 
  group_by(congruency, taskSequence, switchProp, incProp) %>% 
  median_qi(final_drift)
  
#BOUNDARY SEPARATION
bs_qi <- DDM_posteriors %>% 
  group_by(congruency, taskSequence, switchProp, incProp) %>% 
  median_qi(final_bs)

#BOUNDARY SEPARATION
ndt_qi <- DDM_posteriors %>% 
  group_by(congruency, taskSequence, switchProp, incProp) %>% 
  median_qi(final_ndt)

#--------------------------------------------------
# Plots BY CONDITION
drift_posterior_plot <- posterior_plot_func(drift_qi, "final_drift", c("c"="#AA9355", "i"="#843C0C")) +
  theme(legend.position = "none",
        strip.text.y = element_blank())

bs_posterior_plot <- posterior_plot_func(bs_qi, "final_bs", c("c"="#AA9355", "i"="#843C0C")) +
  theme(legend.position = "none",
        strip.text.y = element_blank()) +
  scale_y_continuous(breaks = c(1.5, 1.75, 2))

ndt_posterior_plot <- posterior_plot_func(ndt_qi, "final_ndt", c("c"="#AA9355", "i"="#843C0C")) +
  scale_y_continuous(breaks = c(0.45, .5, .55, .6, 0.65), limits=c(0.44, 0.66))
legend <- cowplot::get_legend(ndt_posterior_plot)
ndt_posterior_plot <- ndt_posterior_plot + 
  theme(legend.position = "none")

#--------------------------------------------------
# get effects for each parameter

#---------------
# DRIFT
drift_effects <- DDM_posteriors %>% 
  group_by(.iteration, .chain, .draw) %>% 
  nest() %>% 
  mutate(coefs = map(data, function(d){coef(lm(final_drift ~ congruency * taskSequence * switchProp * incProp, data=d))})) %>% 
  unnest_wider(coefs) %>% 
  select(-data) %>% 
  rename(!!!rename_list) %>% 
  pivot_longer(cols= -c(.draw, .iteration, .chain),
               names_to = "effect",
               values_to = ".value")

# subselect for this plot, then get qis
drift_plot_effects <- drift_effects %>% 
  filter(effect %in% include_in_plot_list) %>% 
  group_by(effect) %>% 
  median_qi(.value) %>% 
  mutate(sig = case_when(
    .upper > 0 & .lower > 0 ~ 1,
    .upper < 0 & .lower < 0 ~ 1,
    .default = 0
  ))

#fix order
drift_plot_effects$effect = factor(drift_plot_effects$effect, levels = rev(include_in_plot_list))

#---------------
# BS
bs_effects <- DDM_posteriors %>% 
  group_by(.iteration, .chain, .draw) %>% 
  nest() %>% 
  mutate(coefs = map(data, function(d){coef(lm(final_bs ~ congruency * taskSequence * switchProp * incProp, data=d))})) %>% 
  unnest_wider(coefs) %>% 
  select(-data) %>% 
  rename(!!!rename_list) %>% 
  pivot_longer(cols= -c(.draw, .iteration, .chain),
               names_to = "effect",
               values_to = ".value")

# subselect for this plot, then get qis
bs_plot_effects <- bs_effects %>% 
  filter(effect %in% include_in_plot_list) %>% 
  group_by(effect) %>% 
  median_qi(.value) %>% 
  mutate(sig = case_when(
    .upper > 0 & .lower > 0 ~ 1,
    .upper < 0 & .lower < 0 ~ 1,
    .default = 0
  ))

#fix order
bs_plot_effects$effect = factor(bs_plot_effects$effect, levels = rev(include_in_plot_list))

#---------------
# NDT
ndt_effects <- DDM_posteriors %>% 
  group_by(.iteration, .chain, .draw) %>% 
  nest() %>% 
  mutate(coefs = map(data, function(d){coef(lm(final_ndt ~ congruency * taskSequence * switchProp * incProp, data=d))})) %>% 
  unnest_wider(coefs) %>% 
  select(-data) %>% 
  rename(!!!rename_list) %>% 
  pivot_longer(cols= -c(.draw, .iteration, .chain),
               names_to = "effect",
               values_to = ".value")

# subselect for this plot, then get qis
ndt_plot_effects <- ndt_effects %>% 
  filter(effect %in% include_in_plot_list) %>% 
  group_by(effect) %>% 
  median_qi(.value) %>% 
  mutate(sig = case_when(
    .upper > 0 & .lower > 0 ~ 1,
    .upper < 0 & .lower < 0 ~ 1,
    .default = 0
  ))

#fix order
ndt_plot_effects$effect = factor(ndt_plot_effects$effect, levels = rev(include_in_plot_list))

#--------------------------------------------------
# PLOT EFFECTS
drift_effects_plot <- effects_plot_func(drift_plot_effects) 
bs_effects_plot <- effects_plot_func(bs_plot_effects) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ndt_effects_plot <- effects_plot_func(ndt_plot_effects) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
drift_effects_plot + bs_effects_plot + ndt_effects_plot

#--------------------------------------------------
# combine plots

drift_plot <- wrap_elements(drift_posterior_plot / drift_effects_plot)
bs_plot <- wrap_elements(bs_posterior_plot / bs_effects_plot)
ndt_plot <- wrap_elements(ndt_posterior_plot / ndt_effects_plot)

ggsave("plots/drift_pars_replication.jpg", drift_plot, 
       dpi = 300, width = 3.7, height =6, units = "in")
ggsave("plots/bs_pars_replication.jpg", bs_plot, 
       dpi = 300, width = 3, height =6, units = "in")
ggsave("plots/ndt_pars_replication.jpg", ndt_plot, 
       dpi = 300, width = 3.25, height =6, units = "in")

#plot version of drift for y axis labels
ggsave("plots/DDM_weights_legend_replication.jpg", legend,
  dpi = 300, width = 3, height =0.5, units = "in")


