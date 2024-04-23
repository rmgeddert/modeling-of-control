library(tidyverse)
library(tidybayes)
library(plotrix)
library(lme4)
library(kableExtra)
library(stringr)

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
  "Congruency:TaskSequence:SwitchProp:IncProp" = "congruencyi:taskSequences:switchProp75%:incProp75%"
)

include_in_plot_list = c(
  "TaskSequence",
  "Congruency",
  "SwitchProp",
  "IncProp",
  "Congruency:TaskSequence",
  "TaskSequence:SwitchProp",
  "TaskSequence:IncProp",
  "Congruency:SwitchProp",
  "Congruency:IncProp",
  "SwitchProp:IncProp",
  "Congruency:TaskSequence:SwitchProp",
  "Congruency:TaskSequence:IncProp",
  "Congruency:SwitchProp:IncProp",
  "TaskSequence:SwitchProp:IncProp",
  "Congruency:TaskSequence:SwitchProp:IncProp"
)

# -------------------------------------------------------
# REAL DATA

# load real data
df <- read_csv("data/dataset1.csv", col_types = cols()) %>% 
  rename(congruency = stimCongruency,
         taskSequence = switchType) %>% 
  mutate(RT = RT/1000,
         subject = paste0("sub", match(subject, unique(subject))),
         switchProp = ifelse(blockType == "B" | blockType == "D", "75%", "25%"),
         incProp = ifelse(blockType == "A" | blockType == "B", "75%", "25%")) %>% 
  filter(RT < 1.5 & RT > 0.3) %>% 
  select(trialCount, block, acc, RT, congruency, taskSequence, switchProp, incProp, subject)

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
independent_df <- read_csv("model_fitting/2D_model/posterior_pred.csv") 

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
# no tradeoff

notradeoff_df <- read_csv("model_fitting/NoTradeoff_model/posterior_pred.csv") 

notradeoff_effects <- notradeoff_df %>% 
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
notradeoff_plot_effects <- notradeoff_effects %>% 
  filter(effect %in% include_in_plot_list) %>% 
  group_by(effect) %>% 
  median_qi(.value) %>% 
  mutate(model = "No Tradeoff Model")

#fix order
notradeoff_plot_effects$effect = factor(notradeoff_plot_effects$effect, levels = rev(include_in_plot_list))

# ready to plot! :)

# -------------------------------------------------------
# Plotting

# combine independent into single data frame
combined_dfs <- rbind(independent_plot_effects, notradeoff_plot_effects) 

effects_plot <- ggplot(combined_dfs, aes(y = effect, x = .value, color=model)) +
  geom_pointinterval(aes(xmin=.lower, xmax = .upper), position=position_dodge(width=0.7)) +
  geom_point(data = real_plot_effects) +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_discrete(
    labels = c(
      "Congruency:TaskSequence" = "Congruency x\n TaskSequence",
      "TaskSequence:SwitchProp" = "Task Sequence x\n Switch Prop",
      "TaskSequence:IncProp" = "Task Sequence x\n Inc Prop",
      "Congruency:SwitchProp" = "Congruency x\n Switch Prop",
      "Congruency:IncProp" = "Congruency x\n Inc Prop",
      "SwitchProp:IncProp" = "SwitchProp x \n IncProp",
      "Congruency:TaskSequence:SwitchProp" = "Congruency x TaskSequence \nx SwitchProp",
      "Congruency:TaskSequence:IncProp" = "Congruency x TaskSequence \nx IncProp",
      "Congruency:SwitchProp:IncProp" = "Congruency x SwitchProp \nx IncProp",
      "TaskSequence:SwitchProp:IncProp" = "TaskSequence x SwitchProp \nx IncProp",
      "Congruency:TaskSequence:SwitchProp:IncProp" = "Congruency x Task Sequence \nx SwitchProp x IncProp"
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

# # ggsave
# ggsave("plots/effects_plot_ALL.jpg", effects_plot, dpi = 300, width = 6.5, height = 6, units = "in")
# ggsave("plots/effects_plot_legend.jpg", legend, dpi = 300, width = 6.5, height = 0.5, units = "in")


# ------------------------------------------------
# MAKE LATEX TABLE
# save out latex table of posterior median/95CIs and the real means

posterior_table <- combined_dfs %>% 
  mutate(model = ifelse(model == "2D Model", "2D", "1D")) %>% 
  pivot_wider(names_from = model, values_from = c(.value, .lower, .upper), names_sep = ".") %>% 
  left_join(real_plot_effects, by=c("effect")) %>% 
  mutate(.contains.2D = ifelse(.value > .lower.2D & .value  < .upper.2D, 1, 0),
         .contains.1D = ifelse(.value > .lower.1D & .value  < .upper.1D, 1, 0)) %>% 
  select(effect, .value, 
         .value.2D, .lower.2D, .upper.2D, .contains.2D,
         .value.1D, .lower.1D, .upper.1D, .contains.1D) %>% 
  arrange(match(effect, include_in_plot_list))

# mutate posterior_table to add conditional coloring latex text (with cell_spec)
posterior_table <- posterior_table %>% 
  mutate(across(c(".contains.2D", ".contains.1D"), ~ cell_spec(.x, 
                                                               bold = ifelse(.x == 1, TRUE, FALSE),
                                                               background = ifelse(.x == 1, "#5ec962", "white"),
                                                               "latex"))) 

# convert to data.frame, change column names, and print latex
kbl_dat <- data.frame(posterior_table)
row.names(kbl_dat) <- NULL
colnames(kbl_dat) <- c("", "beta", 
                       "median", ".lower", ".upper", "*",
                       "median", ".lower", ".upper", "*")

#create most of the latex table
latex_text <- kbl_dat %>% 
  kbl(digits=3, "latex", escape=FALSE,
      align = "lccccccccc", booktabs = T) %>% 
  add_header_above(c("Effect" = 1, "Real Data" = 1, "2D Model" = 4, "No Tradeoff Model" = 4),
                   bold=T, font_size=12, line=T) %>% 
  collapse_rows(columns = 1:4, latex_hline="none")

# fix some latex manually
change_latex <- function(x, to_remove, replace_with){
  paste(str_replace_all(x, to_remove, replace_with), collapse = "\n")
}


final_latex <- change_latex(latex_text, to_remove = "\\\\cmidrule\\(l\\{3pt\\}r\\{3pt\\}\\)\\{1-1\\}", replace_with = "") %>% 
  change_latex(., to_remove = "5ec962\\}\\{\\\\textbf\\{1\\}", 
               replace_with = "5ec962\\}\\{\\\\textbf\\{\\\\CheckmarkBold\\}") %>% 
  change_latex(., to_remove = "white\\}\\{0\\}", 
               replace_with = "white\\}\\{\\{\\\\fontfamily\\{phv\\}\\\\selectfont\\\\textbf\\{X\\}\\}\\}")

# print out latex with cat to copy paste into overleaf
cat(final_latex, sep="\n")



