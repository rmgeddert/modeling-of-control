library(cmdstanr)
library(tidyverse)
library(tidybayes)
library(bayesplot)
library(bayestestR)
library(patchwork)
library(kableExtra)

setwd("path/to/modeling-of-control")

effects_plot_func <- function(effect_df, title){
  ggplot(effect_df, aes(y = effect, x = .value, xmin=.lower, xmax = .upper, color=source)) +
    geom_pointinterval(position=position_dodge(width=0.7)) +
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
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none") +
    guides(color = guide_legend(override.aes = list(size=10))) +
    scale_color_manual(values = c("#56B1F7", "#132B43"))
}

# for renaming effects of lm on ddm parameters
rename_list = c(
  "Interceept" = "(Intercept)",
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

# list of effects to plot (in case only want to plot some)
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

# function for calculating lme effects
calculate_effects <- function(DDM_posteriors, param){
  effects_df <- DDM_posteriors %>% 
    group_by(.iteration, .chain, .draw) %>% 
    nest() %>% 
    mutate(coefs = map(data, function(d){coef(lm(!! rlang::sym(param) ~ congruency * taskSequence * switchProp * incProp, data=d))})) %>% 
    unnest_wider(coefs) %>% 
    select(-data) %>% 
    rename(!!!rename_list) %>% 
    pivot_longer(cols= -c(.draw, .iteration, .chain),
                 names_to = "effect",
                 values_to = ".value")
  
  # subselect for this plot, then get qis
  plot_effects <- effects_df %>% 
    filter(effect %in% include_in_plot_list) %>% 
    group_by(effect) %>% 
    median_qi(.value) %>% 
    mutate(sig = case_when(
      .upper > 0 & .lower > 0 ~ 1,
      .upper < 0 & .lower < 0 ~ 1,
      .default = 0
    ))
  
  #fix order
  plot_effects$effect = factor(plot_effects$effect, levels = rev(include_in_plot_list))
  
  return(plot_effects)
}

#--------------------------------------------------
# ORIGINAL
og_DDM_posteriors <- read_csv("model_fitting/2D_model/DDM_posteriors.csv")
og_drift_plot_effects <- calculate_effects(og_DDM_posteriors, "final_drift") %>% 
  mutate(source = "Original")
og_bs_plot_effects <- calculate_effects(og_DDM_posteriors, "final_bs") %>% 
  mutate(source = "Original")
og_ndt_plot_effects <- calculate_effects(og_DDM_posteriors, "final_ndt") %>% 
  mutate(source = "Original")

#--------------------------------------------------
# REPLICATION
rep_DDM_posteriors <- read_csv("model_fitting/2D_model_replication/DDM_posteriors.csv")
rep_drift_plot_effects <- calculate_effects(rep_DDM_posteriors, "final_drift") %>% 
  mutate(source = "Replication")
rep_bs_plot_effects <- calculate_effects(rep_DDM_posteriors, "final_bs") %>% 
  mutate(source = "Replication")
rep_ndt_plot_effects <- calculate_effects(rep_DDM_posteriors, "final_ndt") %>% 
  mutate(source = "Replication")

#--------------------------------------------------
# JOIN TOGETHER
drift_plot_effects <- rbind(og_drift_plot_effects, rep_drift_plot_effects)
bs_plot_effects <- rbind(og_bs_plot_effects, rep_bs_plot_effects)
ndt_plot_effects <- rbind(og_ndt_plot_effects, rep_ndt_plot_effects)

#--------------------------------------------------

# PLOT EFFECTS AGAINST EACH OTHER
drift_effects_plot <- effects_plot_func(drift_plot_effects, "Drift Rate") 
bs_effects_plot <- effects_plot_func(bs_plot_effects, "Boundary Separation") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ndt_effects_plot <- effects_plot_func(ndt_plot_effects, "Non-Decision Time") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
drift_effects_plot + bs_effects_plot + ndt_effects_plot

ggsave("plots/2d_pars_original_vs_replication.jpg", drift_effects_plot + bs_effects_plot + ndt_effects_plot, 
       dpi = 300, width = 8, height =6.5, units = "in")


# ------------------------------------------------
# MAKE LATEX TABLE

#make data base
posterior_table <- rbind(mutate(drift_plot_effects, par = "drift"),
                         mutate(bs_plot_effects, par = "bs"),
                         mutate(ndt_plot_effects, par = "ndt")) %>% 
  pivot_wider(names_from = "par", values_from = c(".value", ".lower", ".upper", "sig"), names_sep = ".") %>% 
  select(effect, 
         .value.drift, .lower.drift, .upper.drift, sig.drift,
         .value.bs, .lower.bs, .upper.bs, sig.bs,
         .value.ndt, .lower.ndt, .upper.ndt, sig.ndt) %>% 
  arrange(match(effect, include_in_plot_list))

# add in latex for cells ahead of time
posterior_table <- posterior_table %>% 
  mutate(across(c("sig.drift", "sig.bs", sig.ndt), ~ cell_spec(.x, 
                                                               bold = ifelse(.x == 1, TRUE, FALSE),
                                                               background = ifelse(.x == 1, "#5ec962", "white"),
                                                               "latex"))) 

# make in data.frame and change column names
kbl_dat <- data.frame(posterior_table)
row.names(kbl_dat) <- NULL
colnames(kbl_dat) <- c(" ",
                       "median", ".lower", ".upper", "*",
                       "median", ".lower", ".upper", "*",
                       "median", ".lower", ".upper", "*")

# create latex
latex_text <- kbl_dat %>% 
  kbl(digits=3, "latex", escape=FALSE,
      align = "lcccccccccccc", booktabs = T) %>% 
  add_header_above(c("Effect" = 1, "Drift Rate" = 4, "Boundary Separation" = 4, "Non-Decision Time" = 4),
                   bold=T, font_size=12) 

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


