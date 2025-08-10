#### This script analysis the one-to-many assignments for Parathaa
## Written by Jacob T. Neairng
### Date: Aug 9th 2025


library(dplyr)
library(stringr)
library(patchwork)
setwd("~/Dropbox_Harvard/hutlab/Jacob/Repos/Hills_Project/Parathaa_Project/Revisions2/output_specific_Aug_7_2025/benchmarks/")


### function to create plots
create_multi_assignment_plot <- function(res_df){
  
  par_df <- res_df %>% mutate("# of Parathaa Assignments"=str_count(Species.parathaa, ";") + 1)
  par_df <- par_df %>% mutate("# of NB - Multi Assignments"=str_count(Species.dada, ";") + 1)
  
  ##anything that is NA in Flag's are now incorrect unassigned so we can just call these incorrect.
  par_df$Flag.x[which(is.na(par_df$Flag.x))] <- FALSE
  par_df$Flag.y[which(is.na(par_df$Flag.y))] <- FALSE
  
  ##these are ones that NB didn't assign
  par_df$`# of NB - Multi Assignments`[which(is.na(par_df$`# of NB - Multi Assignments`))] <- 0
  par_df$`# of Parathaa Assignments`[which(is.na(par_df$`# of Parathaa Assignments`))] <- 0
  
  #set up the various types of assignments
  par_df <- par_df %>% 
    mutate(type=case_when(
      Flag.y & Flag.x ~ "Both Correct",
      Flag.y & !Flag.x ~ "Parathaa Correct",
      !Flag.y & Flag.x ~ "NB - Multi Correct",
      !Flag.y & !Flag.x ~ "Both Incorrect"))
  
  
  
  
  
  parathaa_data <- par_df %>% filter(!is.na(Species.parathaa))
  
  
  
  parathaa_plot <- par_df %>% ggplot(aes(x=`# of Parathaa Assignments`, fill=type)) + geom_histogram() +
    theme_bw() + ylab("Number of sequences") +
    ylim(0,8300)
  
  dada2_plot <- par_df %>% ggplot(aes(x=`# of NB - Multi Assignments`, fill=type)) + geom_histogram() +
    theme_bw() + ylab("Number of sequences") +
    ylim(0, 8300)
  
  ##all sequences
  parathaa_plot + dada2_plot + plot_layout(guides="collect")
  
  ggsave("../../../Revisions2/all_seqs_types.png", width=8, height=6)
  
  
  #only ones that got an assignment for each tool
  parathaa_plot2 <- par_df %>% filter(`# of Parathaa Assignments` > 0) %>%
    ggplot(aes(x=`# of Parathaa Assignments`, fill=type)) + geom_histogram() +
    theme_bw() + ylab("Number of sequences")
  
  
  dada2_plot2 <- par_df %>% filter(`# of NB - Multi Assignments` > 0) %>%
    ggplot(aes(x=`# of NB - Multi Assignments`, fill=type)) + geom_histogram() +
    theme_bw() + ylab("Number of sequences")
  
  
  #just those that got at least one assignment
  parathaa_plot2 + dada2_plot2 + plot_layout(guides="collect")
  
  ggsave("../../../Revisions2/atleast_one_types.png", width=8, height=6)    
  
}


create_parathaa_plot <- function(res_df){
  
  par_df <- res_df %>% mutate("# of Parathaa Assignments"=str_count(Species.parathaa, ";") + 1)
  
  ##anything that is NA in Flag's are now incorrect unassigned so we can just call these incorrect.
  par_df$Flag.y[which(is.na(par_df$Flag.y))] <- FALSE
  par_df$`# of Parathaa Assignments`[which(is.na(par_df$`# of Parathaa Assignments`))] <- 0

  
  ##Flag.y indicates if Parathaa was correct or not.
  
  ##filter out things that didn;t get any assignment.
  parathaa_data <- par_df %>% filter(!is.na(Species.parathaa))

  parathaa_plot <- parathaa_data %>% ggplot(aes(x=`# of Parathaa Assignments`, fill=Flag.y)) + geom_histogram(binwidth=1) +
    theme_bw() + ylab("Number of sequences") + labs(fill="Correct Assignment") +
    ylim(c(0,3700)) +
    coord_cartesian(xlim=c(1,12)) +
    scale_x_continuous(breaks=seq(1,12, by=1)) +
    xlab("Number of ambiguous taxa")
    #scale_x_continuous(breaks = seq(floor(min(parathaa_data$`# of Parathaa Assignments`)), 
     #                               ceiling(max(parathaa_data$`# of Parathaa Assignments`)), by = 1))
  
  return(parathaa_plot)
}


## Load V1V2 data
load("holdout1/Figures/synth_mult_arc/V1V2_full_comparisons_adjust.RData")
V1V2_data <- compare.synth_adjust

## Load V4V5 data

load("holdout1/Figures/synth_mult_arc/V4V5_full_comparisons_adjust.RData")
V4V5_data <- compare.synth_adjust


## Load FL data
load("../../output_sens_Aug_7_2025/benchmarks/holdout1/Figures/synth_mult_arc/FL_full_comparisons_adjust.RData")
FL_data_sens <- compare.synth_adjust


V1V2_par_plot <- create_parathaa_plot(V1V2_data)
V4V5_par_plot <- create_parathaa_plot(V4V5_data)
FL_par_plot <- create_parathaa_plot(FL_data_sens)


(V1V2_par_plot + ggtitle("V1V2 Parathaa Specific")) + (V4V5_par_plot + ggtitle("V4V5 Parathaa Specific")) + 
FL_par_plot + ggtitle("FL Parathaa Sensitive") + plot_layout(guides="collect")


