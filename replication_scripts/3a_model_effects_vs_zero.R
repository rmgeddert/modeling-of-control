library(tidyverse)
library(tidybayes)
library(plotrix)
library(lme4)

setwd("path/to/modeling-of-control")

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
  "Congruency:Task Sequence:SwitchProp:IncProp" = "congruencyi:taskSequences:switchProp75%:incProp75%"
)

include_in_plot_list = c(
  "TaskSequence",
  "Congruency",
  "TaskSequence:SwitchProp",
  "TaskSequence:IncProp",
  "Congruency:SwitchProp",
  "Congruency:IncProp"
)

# -------------------------------------------------------
# REAL DATA

# load real data
df <- read_csv("data/dataset2_replication.csv", col_types = cols()) %>% 
  mutate(RT = RT/1000,
         switchProp = ifelse(blockType == "B" | blockType == "D", "75%", "25%"),
         incProp = ifelse(blockType == "A" | blockType == "B", "75%", "25%"),
         subject = match(subject, unique(subject))) %>% 
  filter(RT < 1.5 & RT > 0.3)

df_means <- df %>% 
  group_by(subject, congruency, taskSequence, switchProp, incProp) %>% 
  summarise(RT = mean(RT))

# run a linear model on df_means to get the true beta effects
real_lm <- lm(RT ~ congruency * taskSequence * switchProp * incProp, data = df_means)
real_lmer <- lmer(RT ~ congruency * taskSequence * switchProp * incProp + (1 | subject), data = df)

# get the coefficients
# real_effects <- tibble(as.data.frame(t(coef(real_lm)))) %>% 
real_effects <- tibble(as.data.frame(t(fixef(real_lmer)))) %>% 
  rename(!!!rename_list) %>% 
  pivot_longer(cols = everything(),
               names_to = "effect",
               values_to = ".value") 

# subselect for this plot
real_plot_effects <- real_effects %>% 
  filter(effect %in% include_in_plot_list) %>% 
  mutate(model = "Real Data")

#fix order
real_plot_effects$effect = factor(real_plot_effects$effect, levels = rev(include_in_plot_list))

# ready to plot! :)

# -------------------------------------------------------
# 2D Model ("independent")
independent_df <- read_csv("model_fitting/2D_model_replication/posterior_pred.csv") 

# for each draw, run a linear model, to convert it to the main effects of each thing of interest
independent_effects <- independent_df %>% 
  group_by(.iteration, .chain, .draw) %>% 
  nest() %>% 
  mutate(coefs = map(data, function(d){coef(lm(RT ~ congruency * taskSequence * switchProp * incProp, data=d))})) %>% 
  unnest_wider(coefs) %>% 
  select(-data) %>% 
  rename(!!!rename_list) %>% 
  pivot_longer(cols= -c(.draw, .iteration, .chain),
               names_to = "effect",
               values_to = ".value")

# subselect for this plot, then get qis
independent_plot_effects <- independent_effects %>% 
  filter(effect %in% include_in_plot_list) %>% 
  group_by(effect) %>% 
  median_qi(.value) %>% 
  mutate(model = "2D Model")

#fix order
independent_plot_effects$effect = factor(independent_plot_effects$effect, levels = rev(include_in_plot_list))

# ready to plot! :)

# -------------------------------------------------------
# 1D Model ("tradeoff")

tradeoff_df <- read_csv("model_fitting/1D_model_replication/posterior_pred.csv") 

tradeoff_effects <- tradeoff_df %>% 
  group_by(.iteration, .chain, .draw) %>% 
  nest() %>% 
  mutate(coefs = map(data, function(d){coef(lm(RT ~ congruency * taskSequence * switchProp * incProp, data=d))})) %>% 
  unnest_wider(coefs) %>% 
  select(-data) %>% 
  rename(!!!rename_list) %>% 
  pivot_longer(cols= -c(.draw, .iteration, .chain),
               names_to = "effect",
               values_to = ".value") 

# subselect for this analysis, then get qis
tradeoff_plot_effects <- tradeoff_effects %>% 
  filter(effect %in% include_in_plot_list) %>% 
  group_by(effect) %>% 
  median_qi(.value) %>% 
  mutate(model = "1D Model")

#fix order
tradeoff_plot_effects$effect = factor(tradeoff_plot_effects$effect, levels = rev(include_in_plot_list))

# ready to plot! :)

# -------------------------------------------------------
# Plotting

# combine independent into single data frame
combined_dfs <- rbind(independent_plot_effects, tradeoff_plot_effects) 

effects_plot <- ggplot(combined_dfs, aes(y = effect, x = .value, color=model)) +
  geom_pointinterval(aes(xmin=.lower, xmax = .upper), position=position_dodge(width=0.7)) +
  geom_point(data = real_plot_effects) +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_discrete(
    labels = c(
      "TaskSequence:SwitchProp" = "Task Sequence x\n Switch Prop",
      "TaskSequence:IncProp" = "Task Sequence x\n Inc Prop",
      "Congruency:SwitchProp" = "Congruency x\n Switch Prop",
      "Congruency:IncProp" = "Congruency x\n Inc Prop"
    )
  ) +
  scale_color_manual(name = "Model",
                     values=c("black", "#639C6E", "#9C6391"),
                     breaks=c('Real Data', '2D Model', '1D Model'),
                     labels=c("Real Data", '2D Model', '1D Model')) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(point_size=4, linetype="blank")))
legend <- cowplot::get_legend(effects_plot)
effects_plot <- effects_plot +
  theme(legend.position = "none")

plot(legend)
effects_plot

#save
# ggsave
ggsave("plots/effects_plot_replication.jpg", effects_plot, dpi = 300, width = 6.5, height = 2.5, units = "in")
ggsave("plots/effects_plot_legend_replication.jpg", legend, dpi = 300, width = 6.5, height = 0.5, units = "in")





