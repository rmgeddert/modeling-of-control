library(tidyverse)
library(plotrix)
library(cowplot)
library(patchwork)
library(bayesplot)
library(tidybayes)
library(kableExtra)

setwd("path/to/modeling-of-control")

# ----------------------------------------------------------------------------
# REAL DATA
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
  summarise(sub_RT = mean(RT)) %>% 
  group_by(congruency, taskSequence, switchProp, incProp) %>% 
  summarise(mean_RT = mean(sub_RT), sem_RT = 2*std.error(sub_RT))

# ----------------------------------------------------------------------------
# Load posteriors and get median/95% CI

posterior_2d <- read_csv("model_fitting/2D_model/posterior_pred.csv") %>% 
  group_by(congruency, taskSequence, switchProp, incProp) %>% 
  median_qi(RT)

posterior_1d <- read_csv("model_fitting/1D_model/posterior_pred.csv") %>%
  group_by(congruency, taskSequence, switchProp, incProp) %>%
  median_qi(RT)

# ----------------------------------------------------------------------------
# Make  Plot

# for dodgning the real data
position_nudgedodge <- function(x = 0, y = 0, width = 0.75) {
  ggproto(NULL, PositionNudgedodge,
          x = x,
          y = y,
          width = width
  )
}

PositionNudgedodge <- ggproto("PositionNudgedodge", PositionDodge,
                              x = 0,
                              y = 0,
                              width = 0.3,
                              
                              setup_params = function(self, data) {
                                l <- ggproto_parent(PositionDodge,self)$setup_params(data)
                                append(l, list(x = self$x, y = self$y))
                              },
                              
                              compute_layer = function(self, data, params, layout) {
                                d <- ggproto_parent(PositionNudge,self)$compute_layer(data,params,layout)
                                d <- ggproto_parent(PositionDodge,self)$compute_layer(d,params,layout)
                                d
                              }
)

#for nice nudging
df_means <- df_means %>%
  mutate(nudge_direction = ifelse(taskSequence == "s", 0.2, -0.2))

plot_func <- function(posterior_df, color_labels){
  pd <- position_dodge(0.2)
  
  ggplot(posterior_df, aes(x=taskSequence, y=RT, color=congruency, group=congruency, shape=congruency)) +
    geom_pointinterval(aes(ymin=.lower, ymax=.upper), position=pd) +
    geom_errorbar(data=df_means, aes(x=taskSequence, y=mean_RT, ymin=mean_RT - sem_RT, ymax=mean_RT + sem_RT),
                  position = position_nudgedodge(x = df_means$nudge_direction, width=0.2),
                  color="black",
                  width=0.1) +
    geom_point(data=df_means, aes(x=taskSequence, y=mean_RT),
               position = position_nudgedodge(x = df_means$nudge_direction, width=0.2),
               color="black", show.legend = FALSE) +
    geom_line(position=pd, show.legend = FALSE) +
    facet_grid(switchProp ~ incProp,
               labeller = labeller(incProp = c("25%" = "25% Incongruent",
                                               "75%" = "75% Incongruent"),
                                   switchProp = c("25%" = "25% Switch",
                                                  "75%" = "75% Switch"))) +
    scale_color_manual(name="Congruency",
                       values = color_labels,
                       labels = c("c"="Congruent", "i"="Incongruent")) +
    scale_shape_manual(name = "Congruency",
                       values = c("c" = 15, "i" = 17),
                       labels = c("c"="Congruent", "i"="Incongruent")) +
    scale_x_discrete(labels=c("r" = "Repeat", "s" = "Switch")) +
    xlab("Task Sequence") +
    guides(color = guide_legend(override.aes = list(point_size=4,
                                                    linetype="blank")))
}

#plot layout:
#556CAA <- nice blue
# yaxis label -> 2d -> 1d 
# _____ -> xaxis label -> xaxis_label
plot_2d <- plot_func(posterior_2d, c("c"="#AA9355", "i"="#843C0C")) +
  theme(legend.position = "none",
        strip.text.y = element_blank()) 
plot_1d <- plot_func(posterior_1d, c("c"="#AA9355", "i"="#843C0C")) +
  theme(legend.position="bottom",
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
legend <- cowplot::get_legend(plot_1d)
plot_1d <- plot_1d + theme(legend.position="none")

plot(legend)
plot_2d + plot_1d

ggsave("plots/posteriors_plot.jpg", plot_2d + plot_1d, dpi = 300, width = 6.5, height = 4, units = "in")
ggsave("plots/posteriors_plot_legend.jpg", legend, dpi = 300, width = 6.5, height = 0.5, units = "in")

# ------------------------------------------------
# MAKE LATEX TABLE
# save out latex table of posterior median/95CIs and the real means
posterior_table <- posterior_2d %>% 
  left_join(posterior_1d, by=setdiff(names(posterior_2d), c("RT", ".lower", ".upper")), suffix=c(".2D", ".1D")) %>% 
  left_join(df_means, by=c("congruency", "taskSequence", "switchProp", "incProp")) %>% 
  mutate(.contains.2D = ifelse(mean_RT > .lower.2D & mean_RT < .upper.2D, 1, 0),
         .contains.1D = ifelse(mean_RT > .lower.1D & mean_RT < .upper.1D, 1, 0)) %>% 
  select(switchProp, incProp, congruency, taskSequence, mean_RT, 
         RT.2D, .lower.2D, .upper.2D, .contains.2D,
         RT.1D, .lower.1D, .upper.1D, .contains.1D) %>%
  arrange(switchProp, incProp) 

# mutate posterior_table to add conditional coloring latex text (with cell_spec)
posterior_table <- posterior_table %>% 
  mutate(across(c(switchProp, incProp), ~ ifelse(.x == "25%", "25\\%", "75\\%"))) %>% 
  mutate(across(c(".contains.2D", ".contains.1D"), ~ cell_spec(.x, 
                                                               bold = ifelse(.x == 1, TRUE, FALSE),
                                                               background = ifelse(.x == 1, "#5ec962", "white"),
                                                               "latex"))) 

# convert to data.frame, change column names, and print latex
kbl_dat <- data.frame(posterior_table)
row.names(kbl_dat) <- NULL
colnames(kbl_dat) <- c("switchProp", "incProp", "congruency", "taskSequence", "mean", 
                                              "median", ".lower", ".upper", "*",
                                              "median", ".lower", ".upper", "*")

latex_text <- kbl_dat %>% 
  kbl(digits=3, "latex", escape=FALSE,
      align = "ccccccccccccc", booktabs = T) %>% 
  add_header_above(c("Conditions" = 4, "Real Data" = 1, "2D Model" = 4, "1D Model" = 4),
                   bold=T, font_size=12) %>%  
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
