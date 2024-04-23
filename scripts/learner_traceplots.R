library(patchwork)
library(Dict)
library(cowplot)

setwd("path/to/modeling-of-control")

blockMapping <- Dict$new(
  "A" = c(0.25, 0.75),
  "B" = c(0.75, 0.75),
  "C" = c(0.25, 0.25),
  "D" = c(0.75, 0.25)
)

#' Creates a task data set for a single participant, drawing on switch and incongruent proportion
#' for each block based on given blockTypeOrder and the blockMapping
#' 
#' @param nBlocks blocks in the task, should be 4
#' @param nTrialsPerBlock trials in each block, should be 128
#' @param blockTypeOrder Order of blocks, usually counterbalanced with latin square
createTaskData <- function(nBlocks=4,nTrialsPerBlock=128,blockTypeOrder=c("A", "B", "D", "C")){
  # vector of each trials task sequence, based on block probabilities and given block order
  taskSequence = c(
    "n", sample(c(rep("s",nTrialsPerBlock * blockMapping[blockTypeOrder[1]][1] - 1), 
                  rep("r",nTrialsPerBlock * (1 - blockMapping[blockTypeOrder[1]][1])))),
    "n", sample(c(rep("s",nTrialsPerBlock * blockMapping[blockTypeOrder[2]][1] - 1), 
                  rep("r",nTrialsPerBlock * (1 - blockMapping[blockTypeOrder[2]][1])))),
    "n", sample(c(rep("s",nTrialsPerBlock * blockMapping[blockTypeOrder[3]][1] - 1), 
                  rep("r",nTrialsPerBlock * (1 - blockMapping[blockTypeOrder[3]][1])))),
    "n", sample(c(rep("s",nTrialsPerBlock * blockMapping[blockTypeOrder[4]][1] - 1), 
                  rep("r",nTrialsPerBlock * (1 - blockMapping[blockTypeOrder[4]][1]))))
  )
  
  # vector of each trials congruency, based on block probabilities
  congruency = c(
    sample(c(rep("i",nTrialsPerBlock * blockMapping[blockTypeOrder[1]][2]), 
             rep("c",nTrialsPerBlock * (1 - blockMapping[blockTypeOrder[1]][2])))),
    sample(c(rep("i",nTrialsPerBlock * blockMapping[blockTypeOrder[2]][2]), 
             rep("c",nTrialsPerBlock * (1 - blockMapping[blockTypeOrder[2]][2])))),
    sample(c(rep("i",nTrialsPerBlock * blockMapping[blockTypeOrder[3]][2]), 
             rep("c",nTrialsPerBlock * (1 - blockMapping[blockTypeOrder[3]][2])))),
    sample(c(rep("i",nTrialsPerBlock * blockMapping[blockTypeOrder[4]][2]), 
             rep("c",nTrialsPerBlock * (1 - blockMapping[blockTypeOrder[4]][2]))))
  )
  
  # create tibble of task
  task_data <- tibble(
    trialCount = 1:(nTrialsPerBlock*nBlocks),
    block = rep(1:nBlocks, each=nTrialsPerBlock),
    blockType = blockTypeOrder[block],
    switchProp = ifelse(blockType == "A" | blockType == "C", 0.25, 0.75),
    incProp =ifelse(blockType == "A" | blockType == "B", 0.75, 0.25),
    taskSequence = taskSequence,
    congruency = congruency,
    ts_val = ifelse(taskSequence == "n", 0.5, ifelse(taskSequence == "s", 1, 0)),
    inc_val = ifelse(congruency == "i", 1, 0),
  )
}

# create fake task data
df <- createTaskData() 

# apply ideal proportion learner
df <- df %>% 
  group_by(block) %>% 
  mutate(
    inc_sum = lag(cumsum(inc_val), default=0),
    switch_sum = lag(cumsum(ts_val), default=0),
    count = row_number(),
    switchPropApprox = ifelse(count == 1, 0.5, switch_sum / (count - 1)),
    incPropApprox = ifelse(count == 1, 0.5, inc_sum / (count - 1))
  )

# add stability, flexibility, and stability-flexibility level
df <- df %>% 
  mutate(
    stability=incPropApprox,
    flexibility = switchPropApprox,
    stability_flexibility = (incPropApprox + (1 - switchPropApprox)) / 2
  )

p1 <- ggplot(df, aes(x=trialCount)) +
  geom_line(aes(y=switchPropApprox, color="Switch"), size=0.8) + 
  geom_line(aes(y=incPropApprox, color="Incongruency"), size=0.8) +
  geom_vline(xintercept = seq(128, 384, by=128)) +
  scale_y_continuous(labels=c("0","0.25","0.50","0.75","1.00"),
                     breaks=c(0,0.25, 0.5, 0.75, 1),
                     expand = c(0,0)) +
  scale_x_continuous(limits = c(0, 512),
                     expand = c(0,0),
                     breaks = c(0, 100, 200, 300, 400, 500),
                     labels=c("0", "100", "200", "300", "400", "500")) +
  labs(y = "Proportion") +
  scale_color_manual(values = c("Switch" = "#468D30", "Incongruency" = "#FF8D30")) +
  ggtitle("Proportion Estimates") +
  theme_cowplot() +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.justification = "center",
        plot.title = element_text(hjust = 0.5))

p2 <- ggplot(df, aes(x=trialCount)) +
  geom_line(aes(y=stability, color="Stability"), size=0.8) + 
  geom_line(aes(y=flexibility, color="Flexibility"), size=0.8) +
  geom_vline(xintercept = seq(128, 384, by=128)) +
  scale_y_continuous(labels=c("0","0.25","0.50","0.75","1.00"),
                     breaks=c(0,0.25, 0.5, 0.75, 1),
                     expand = c(0,0)) +
  scale_x_continuous(limits = c(0, 512),
                     expand = c(0,0),
                     breaks = c(0, 100, 200, 300, 400, 500),
                     labels=c("0", "100", "200", "300", "400", "500")) +
  expand_limits(y=c(0,1)) +
  labs(x = "Trial Count") +
  scale_color_manual(values = c("Stability" = "#FF8D30", "Flexibility" = "#468D30"),
                     breaks = c("Stability", "Flexibility")) +
  ggtitle("2D Model") +
  theme_cowplot() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.justification = "center",
        plot.title = element_text(hjust = 0.5)) 

p3 <- ggplot(df, aes(x=trialCount, y=stability_flexibility)) +
  geom_line(aes(color = "Stability-Flexibility"), size=0.8) + 
  geom_vline(xintercept = seq(128, 384, by=128)) +
  scale_y_continuous(labels=c("0","0.25","0.50","0.75","1.00"),
                     breaks=c(0,0.25, 0.5, 0.75, 1),
                     expand = c(0,0)) +
  scale_x_continuous(limits = c(0, 512),
                     expand = c(0,0),
                     breaks = c(0, 100, 200, 300, 400, 500),
                     labels=c("0", "100", "200", "300", "400", "500")) +
  expand_limits(y=c(0,1)) +
  labs(x = "Trial Count") +
  scale_color_manual(labels = c("Stability-Flexibility"),
                     values = c("Stability-Flexibility" = "#7BD3EA")) +
  ggtitle("1D Model") +
  theme_cowplot() +
  theme(axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.justification = "center",
        plot.title = element_text(hjust = 0.5))

p1 / p2 / p3 

ggsave("plots/learner_traceplots.jpg", p1 / p2 / p3, dpi = 300, width = 8, height = 10, units = "in")
