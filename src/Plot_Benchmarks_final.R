#### Plot benchmark figures
require(docopt)

'Usage:
  Plot_Benchmarks.R [--specific <benchmarks_folder_for_specifc> --sensitive <benchmark_data_for_sensitive> --output <output_directory>]
  
  
Options:
  --specific The folder with the specific (default) benchmark results
  --sensitive The folder with the sensitive benchmark results
  --output the output directory for plots
  ]' -> doc


opts <- docopt(doc)

opts$specific <- "~/Dropbox_Harvard/hutlab/Jacob/Repos/Hills_Project/Parathaa_Project/Main_results_July_22/output_specific/benchmarks/"
opts$sensitive <- "~/Dropbox_Harvard/hutlab/Jacob/Repos/Hills_Project/Parathaa_Project/Main_results_July_22/output/benchmarks/"
opts$output <- "~/Desktop"
library(dplyr)
library(ggplot2)
library(reshape2)

## function that takes in the folder and generates list of dataframes for each benchmark and region
get_data <- function(root_dir, level, adjust_name="", mult=""){
  
  V1V2_res <- list()
  V4V5_res <- list()
  Full_res_naive<- list()
  Full_res_exact <- list()
  
  directories <- list.dirs(root_dir, recursive = F)
  
  
  V1V2_extension <- paste0("/Figures/synth_", mult, "_arc/V1V2_", level, "_performance_", adjust_name, ".tsv")
  V4V5_extension <- paste0("/Figures/synth_", mult, "_arc/V4V5_", level, "_performance_", adjust_name, ".tsv")
  Full_extension_naive <- paste0("/FL/FL_", level, "_performance_", adjust_name, ".tsv")
  Full_extension_exact <- paste0("/Figures/synth_", mult, "_arc/FL_", level, "_performance_", adjust_name, ".tsv")
  
  for(dir in directories){
    name <- gsub(".*\\/", "", dir)
    V1V2_res[[paste0(name, mult)]] <- read.table(paste0(dir, V1V2_extension))
    V4V5_res[[paste0(name, mult)]] <- read.table(paste0(dir, V4V5_extension))
    Full_res_naive[[name]] <- read.table(paste0(dir, Full_extension_naive))
    Full_res_exact[[name]] <- read.table(paste0(dir, Full_extension_exact))
  }
  
  final_results <- list(V1V2=V1V2_res, V4V5=V4V5_res, Full_naive=Full_res_naive, Full_exact=Full_res_exact)
  
  return(final_results)  
  
}

#converts the list of dataframes into a single melted dataframe
melt_bench_data <- function(bench_data){
  region_dfs <- list()
  
  for(i in 1:length(bench_data)){
    region_dfs[[names(bench_data)[i]]] <- do.call(rbind, bench_data[[i]])
  }
  
  full_df <- do.call(rbind, region_dfs)
  full_df$dataset <- sub("\\.[0-9]+", "", rownames(full_df))
  full_df$dataset <- sub(".*\\.", "", full_df$dataset)
  full_df$dataset <- sub("mult", "", full_df$dataset)
  full_df$region <- sub("\\..*", "", rownames(full_df))
  
  colnames(full_df) <- c("Metric", "Parathaa", "DADA2", "Dataset", "Region")
  
  #remove the rows with blank metric
  rm_index <- which(full_df$Metric=="")
  
  full_df <- full_df[-rm_index,]
  
  full_df_melt <- melt(full_df, id.vars = c("Metric", "Dataset", "Region"))
  
  
  return(full_df_melt)
}

plot_benchs_species <- function(mult_bench, no_mult_bench, sens_bench){
  
  #load in each data frame 
  mult_df_melt <- melt_bench_data(mult_bench)
  no_mult_df_melt <- melt_bench_data(no_mult_bench)
  sens_bench_df_melt <- melt_bench_data(sens_bench)
  
  
  #Add Parathaa specific results from Full, V1V2
  no_mult_df_melt$Dataset <- gsub("no$", "", no_mult_df_melt$Dataset)
  #remove full naive Parathaa since its the same result as full exact (i.e. nothing changes between them)
  no_mult_df_melt <- no_mult_df_melt %>% filter(Region != "Full_naive" | variable != "Parathaa")
  no_mult_df_melt$variable <- as.character(no_mult_df_melt$variable)
  no_mult_df_melt$variable[which(no_mult_df_melt$Region=="Full_naive" & no_mult_df_melt$variable=="DADA2")] <- "DADA2 Naive Bayes"
  no_mult_df_melt$variable[which(no_mult_df_melt$Region=="Full_exact" & no_mult_df_melt$variable=="DADA2")] <- "DADA2 Exact Match"
  no_mult_df_melt$Region <- gsub("_.*", "", no_mult_df_melt$Region)
  
  
  
  
  
  #add dada2 multi results from V1V2 and V4V5 regions
  dada2_mult <- mult_df_melt %>% filter(variable=="DADA2" & Region != "Full_exact" & Region != "Full_naive")
  dada2_mult$variable <- "DADA2 Multi"  
  full_df_melt <- rbind(no_mult_df_melt, dada2_mult)
  
  #add sensitive parathaa results
  sens_bench_df_melt <- sens_bench_df_melt %>% filter(variable=="Parathaa") %>% 
    filter(Region %in% c("V1V2", "V4V5", "Full_naive"))
  sens_bench_df_melt$variable <- "Parathaa Sensitive"
  sens_bench_df_melt$Region <- gsub("_.*", "", sens_bench_df_melt$Region)
  
  full_df_melt <- rbind(full_df_melt, sens_bench_df_melt)
  full_df_melt$variable <- as.character(full_df_melt$variable)
  
  full_df_melt$variable[which(full_df_melt$variable=="Parathaa")] <- "Parathaa Specific (default)"
  
  full_df_melt$Metric <- factor(full_df_melt$Metric, levels=c("F1 Score", "Accuracy", "Precision", 
                                                              "Recall", "Uniquely Correct",
                                                              "One-to-many Correct", "Incorrect", "Unassigned Correct", "Unassigned Incorrect"))
  
  
  
  cols <- c("Parathaa Specific (default)"="#FA8072",
            "Parathaa Sensitive"="#660000",
            "DADA2"="#143d59",
            "DADA2 Multi"="#1e1eff",
            "DADA2 Naive Bayes"="#47008e",
            "DADA2 Exact Match"="#82EEFD")
  
  shapes <- c("Parathaa Specific (default)"=16,
              "Parathaa Sensitive"=16,
              "DADA2"=17,
              "DADA2 Multi"=17,
              "DADA2 Naive Bayes"=17,
              "DADA2 Exact Match"=17)
  
  full_df_melt$Dataset <- gsub("holdout[1-3]", "Holdout", full_df_melt$Dataset)
  full_df_melt$Dataset <- gsub("original", "Baseline", full_df_melt$Dataset)
  full_df_melt$Dataset <- gsub("holdoutOG", "Baseline Holdout", full_df_melt$Dataset)
  full_df_melt$Dataset <- gsub("novel", "Genus outside of DB", full_df_melt$Dataset)
  full_df_melt$Dataset <- gsub("even", "Single represenative genus", full_df_melt$Dataset)
  
  full_df_melt$variable <- factor(full_df_melt$variable, levels=c("Parathaa Specific (default)", 
                                                                  "DADA2 Multi",
                                                                  "DADA2 Exact Match",
                                                                  "Parathaa Sensitive",
                                                                  "DADA2",
                                                                  "DADA2 Naive Bayes"
                                                                  ))
  
  full_df_melt$Dataset <- factor(full_df_melt$Dataset, levels=c("Baseline", "Baseline Holdout", "Holdout",
                                                                "Single represenative genus",  "Genus outside of DB"))
  
  full_plot <- ggplot(full_df_melt, aes(x=Dataset, y=as.numeric(value), group=variable, color=variable, shape=variable)) + 
    geom_jitter(position=position_jitterdodge(dodge.width = 1, jitter.width = .1), size=3, alpha=0.7) +
    #geom_point(stat="summary", fun="mean", shape="_", position = position_dodge(width = 1), size=10) +
    facet_grid(rows=vars(Metric), cols=vars(Region), scales="free_y") +
    theme_bw() + 
    scale_color_manual(values = cols, name="Software", breaks=
                         c("Parathaa Specific (default)", "Parathaa Sensitive", "DADA2", "DADA2 Multi", 
                           "DADA2 Exact Match", "DADA2 Naive Bayes"))+ 
    ylab("") +
    ggtitle("Species") + theme(strip.text.y=element_text(angle=0), axis.text.x = element_text(angle=65, vjust = 0.5)) +
    scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
    scale_shape_manual(name="Software", values=shapes, breaks=
                         c("Parathaa Specific (default)", "Parathaa Sensitive", "DADA2", "DADA2 Multi", 
                           "DADA2 Exact Match", "DADA2 Naive Bayes")) +
    geom_rect(xmin=0,
              xmax=1.5,
              ymin=0, ymax=1, color=NA, fill="grey", alpha=0.01) +
    geom_rect(xmin=1.5,
              xmax=2.5,
              ymin=0, ymax=1, color=NA, fill="black", alpha=0.01) +
    geom_rect(xmin=2.5,
              xmax=3.5,
              ymin=0, ymax=1, color=NA, fill="grey", alpha=0.01) +
    geom_rect(xmin=3.5,
              xmax=4.5,
              ymin=0, ymax=1, color=NA, fill="black", alpha=0.01) +
    geom_rect(xmin=4.5,
              xmax=6,
              ymin=0, ymax=1, color=NA, fill="grey", alpha=0.01)
  
  filt_data <- full_df_melt %>% filter(Dataset=="Holdout") %>% filter(Metric %in% c("F1 Score", "Precision", "Recall", 
                                                                                    "Uniquely Correct", "One-to-many Correct",
                                                                                    "Incorrect"))
  
  top_panel <- c("F1 Score", "Precision", "Recall")
  bottom_panel <- c("Uniquely Correct", "One-to-many Correct",
                    "Incorrect")
  dumby_data <- filt_data
  dumby_data <- dumby_data %>% mutate(value = case_when(
    Metric %in% top_panel & variable == "Parathaa Specific (default)" ~ 1,
    Metric %in% top_panel & variable == "DADA2" ~ 0.4,
    Metric %in% top_panel & variable == "DADA2 Exact Match" ~ 0.4,
    Metric %in% bottom_panel & variable == "Parathaa Specific (default)" ~ 0.5,
    Metric %in% bottom_panel & variable == "DADA2" ~ 0.5,
    Metric %in% bottom_panel & variable == "DADA2 Naive Bayes" ~ 0.5,
  ))
  
  small_plot <- ggplot(filt_data, aes(x=Dataset, y=as.numeric(value), group=variable, color=variable, shape=variable)) + 
    geom_jitter(position=position_jitterdodge(dodge.width = 1, jitter.width = .2), size=3, alpha=0.6) +
    #geom_point(stat="summary", fun="mean", shape="_", position = position_dodge(width = 1), size=10) +
    facet_grid(rows=vars(Metric), cols=vars(Region), scales="free_y") +
    theme_bw() + 
    scale_color_manual(values = cols, name="Software", breaks=
                         c("Parathaa Specific (default)", "Parathaa Sensitive", "DADA2", "DADA2 Multi", 
                           "DADA2 Exact Match", "DADA2 Naive Bayes"))+ 
    ylab("") +
    ggtitle("Species") + theme(strip.text.y=element_text(angle=0), axis.text.x = element_text(angle=65, vjust = 0.5)) +
    scale_shape_manual(values=shapes, name="Software", breaks=
                         c("Parathaa Specific (default)", "Parathaa Sensitive", "DADA2", "DADA2 Multi", 
                           "DADA2 Exact Match", "DADA2 Naive Bayes"))+
    scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
    geom_blank(data=dumby_data) +
    theme(panel.spacing.y = unit(.5, "lines")) +
    xlab("") +
    theme(axis.text.x=element_blank(), axis.ticks.x =  element_blank())
  
  
  return(list(full_plot, small_plot))
  
}


mult_adjusted_Species <- get_data(opts$specific, level="Species", adjust_name="adjust", mult="mult")
nomult_adjusted_Species <- get_data(opts$specific, level="Species", adjust_name = "adjust", mult="nomult")
sensitive_Species <- get_data(opts$sensitive, level="Species", adjust_name = "adjust", mult="mult")




species_plots <- plot_benchs_species(mult_bench = mult_adjusted_Species,
                                     no_mult_bench = nomult_adjusted_Species,
                                     sens_bench = sensitive_Species)




library(grid)
gt = ggplot_gtable(ggplot_build(species_plots[[2]]))
gt$heights[15] = 3*gt$heights[15]
grid.draw(gt)

ggsave(filename = file.path(opts$output, "Species_full_data.pdf"), plot=species_plots[[1]], height = 8.25, width=11, units = "in")
ggsave(filename = file.path(opts$output, "Species_filt_data.pdf"), plot=species_plots[[2]], height = 6, width=8, units = "in")


##### Revisions


species_plots[[2]]$data$variable <- as.character(species_plots[[2]]$data$variable)

## change names

species_plots[[2]]$data$variable[which(species_plots[[2]]$data$variable=="DADA2 Multi")] <- "Naive Bayes Multi-EM"
species_plots[[2]]$data$variable[which(species_plots[[2]]$data$variable=="DADA2 Exact Match")] <- "Naive Bayes Multi-EM"
species_plots[[2]]$data$variable[which(species_plots[[2]]$data$variable=="DADA2 Naive Bayes")] <- "Naive Bayes"
species_plots[[2]]$data$variable[which(species_plots[[2]]$data$variable=="DADA2")] <- "Naive Bayes-EM"


species_plots[[1]]$data$variable <- as.character(species_plots[[1]]$data$variable)
species_plots[[1]]$data$variable[which(species_plots[[1]]$data$variable=="DADA2 Multi")] <- "Naive Bayes Multi-EM"
species_plots[[1]]$data$variable[which(species_plots[[1]]$data$variable=="DADA2 Exact Match")] <- "Naive Bayes Multi-EM"
species_plots[[1]]$data$variable[which(species_plots[[1]]$data$variable=="DADA2 Naive Bayes")] <- "Naive Bayes"
species_plots[[1]]$data$variable[which(species_plots[[1]]$data$variable=="DADA2")] <- "Naive Bayes-EM"


plot_data <- species_plots[[2]]$data


### add in TaxaID data
FL_id_taxa <- read.table("~/Dropbox_Harvard/hutlab/Jacob/Repos/Hills_Project/Parathaa_Project/IDTaxa_testing/FINAL_RESULTS/holdout1/FL_Species_performance_adjust.tsv", header=T)

FL_id_taxa <- FL_id_taxa %>% melt()
colnames(FL_id_taxa)[1] <- c("Metric")
FL_id_taxa$Region <- "Full"
FL_id_taxa$variable <- "IDTAXA"
colnames(FL_id_taxa)
FL_id_taxa$Dataset <- "Holdout"
FL_id_taxa1 <- FL_id_taxa %>% filter(!Metric %in% c("Unassigned Correct", "Unassigned Incorrect", "Accuracy"))
FL_id_taxa1_full <- FL_id_taxa

FL_id_taxa <- read.table("~/Dropbox_Harvard/hutlab/Jacob/Repos/Hills_Project/Parathaa_Project/IDTaxa_testing/FINAL_RESULTS/IDTaxa/holdout2/FL_Species_performance_adjust.tsv", header=T)

FL_id_taxa <- FL_id_taxa %>% melt()
colnames(FL_id_taxa)[1] <- c("Metric")
FL_id_taxa$Region <- "Full"
FL_id_taxa$variable <- "IDTAXA"
colnames(FL_id_taxa)
FL_id_taxa$Dataset <- "Holdout"
FL_id_taxa2 <- FL_id_taxa %>% filter(!Metric %in% c("Unassigned Correct", "Unassigned Incorrect", "Accuracy"))
FL_id_taxa2_full <- FL_id_taxa

FL_id_taxa <- read.table("~/Dropbox_Harvard/hutlab/Jacob/Repos/Hills_Project/Parathaa_Project/IDTaxa_testing/FINAL_RESULTS/IDTaxa/holdout3//FL_Species_performance_adjust.tsv", header=T)

FL_id_taxa <- FL_id_taxa %>% melt()
colnames(FL_id_taxa)[1] <- c("Metric")
FL_id_taxa$Region <- "Full"
FL_id_taxa$variable <- "IDTAXA"
colnames(FL_id_taxa)
FL_id_taxa$Dataset <- "Holdout"
FL_id_taxa3 <- FL_id_taxa %>% filter(!Metric %in% c("Unassigned Correct", "Unassigned Incorrect", "Accuracy"))
FL_id_taxa3_full <- FL_id_taxa

FL_id_taxa <- rbind(FL_id_taxa1, FL_id_taxa2, FL_id_taxa3)

## load in V1V2
V1V2_id_taxa <- read.table("~/Dropbox_Harvard/hutlab/Jacob/Repos/Hills_Project/Parathaa_Project/IDTaxa_testing/FINAL_RESULTS/holdout1/V1V2_Species_performance_adjust.tsv", header=T)

V1V2_id_taxa <- V1V2_id_taxa %>% melt()
colnames(V1V2_id_taxa)[1] <- c("Metric")
V1V2_id_taxa$Region <- "V1V2"
V1V2_id_taxa$variable <- "IDTAXA"
colnames(V1V2_id_taxa)
V1V2_id_taxa$Dataset <- "Holdout"
V1V2_id_taxa1 <- V1V2_id_taxa %>% filter(!Metric %in% c("Unassigned Correct", "Unassigned Incorrect", "Accuracy"))
V1V2_id_taxa1_full <- V1V2_id_taxa


V1V2_id_taxa <- read.table("~/Dropbox_Harvard/hutlab/Jacob/Repos/Hills_Project/Parathaa_Project/IDTaxa_testing/FINAL_RESULTS/IDTaxa/holdout2//V1V2_Species_performance_adjust.tsv", header=T)

V1V2_id_taxa <- V1V2_id_taxa %>% melt()
colnames(V1V2_id_taxa)[1] <- c("Metric")
V1V2_id_taxa$Region <- "V1V2"
V1V2_id_taxa$variable <- "IDTAXA"
colnames(V1V2_id_taxa)
V1V2_id_taxa$Dataset <- "Holdout"
V1V2_id_taxa2 <- V1V2_id_taxa %>% filter(!Metric %in% c("Unassigned Correct", "Unassigned Incorrect", "Accuracy"))
V1V2_id_taxa2_full <- V1V2_id_taxa

V1V2_id_taxa <- read.table("~/Dropbox_Harvard/hutlab/Jacob/Repos/Hills_Project/Parathaa_Project/IDTaxa_testing/FINAL_RESULTS/IDTaxa/holdout3//V1V2_Species_performance_adjust.tsv", header=T)

V1V2_id_taxa <- V1V2_id_taxa %>% melt()
colnames(V1V2_id_taxa)[1] <- c("Metric")
V1V2_id_taxa$Region <- "V1V2"
V1V2_id_taxa$variable <- "IDTAXA"
colnames(V1V2_id_taxa)
V1V2_id_taxa$Dataset <- "Holdout"
V1V2_id_taxa3 <- V1V2_id_taxa %>% filter(!Metric %in% c("Unassigned Correct", "Unassigned Incorrect", "Accuracy"))
V1V2_id_taxa3_full <- V1V2_id_taxa

V1V2_id_taxa <- rbind(V1V2_id_taxa1, V1V2_id_taxa2, V1V2_id_taxa3)

## Load in V4V5

V4V5_id_taxa <- read.table("~/Dropbox_Harvard/hutlab/Jacob/Repos/Hills_Project/Parathaa_Project/IDTaxa_testing/FINAL_RESULTS/holdout1/V4V5_Species_performance_adjust.tsv", header=T)
V4V5_id_taxa <- V4V5_id_taxa %>% melt()
colnames(V4V5_id_taxa)[1] <- c("Metric")
V4V5_id_taxa$Region <- "V4V5"
V4V5_id_taxa$variable <- "IDTAXA"
colnames(V4V5_id_taxa)
V4V5_id_taxa$Dataset <- "Holdout"
V4V5_id_taxa1 <- V4V5_id_taxa %>% filter(!Metric %in% c("Unassigned Correct", "Unassigned Incorrect", "Accuracy"))
V4V5_id_taxa1_full <- V4V5_id_taxa

V4V5_id_taxa <- read.table("~/Dropbox_Harvard/hutlab/Jacob/Repos/Hills_Project/Parathaa_Project/IDTaxa_testing/FINAL_RESULTS/IDTaxa/holdout2//V4V5_Species_performance_adjust.tsv", header=T)
V4V5_id_taxa <- V4V5_id_taxa %>% melt()
colnames(V4V5_id_taxa)[1] <- c("Metric")
V4V5_id_taxa$Region <- "V4V5"
V4V5_id_taxa$variable <- "IDTAXA"
colnames(V4V5_id_taxa)
V4V5_id_taxa$Dataset <- "Holdout"
V4V5_id_taxa2 <- V4V5_id_taxa %>% filter(!Metric %in% c("Unassigned Correct", "Unassigned Incorrect", "Accuracy"))
V4V5_id_taxa2_full <- V4V5_id_taxa

V4V5_id_taxa <- read.table("~/Dropbox_Harvard/hutlab/Jacob/Repos/Hills_Project/Parathaa_Project/IDTaxa_testing/FINAL_RESULTS/IDTaxa/holdout3//V4V5_Species_performance_adjust.tsv", header=T)
V4V5_id_taxa <- V4V5_id_taxa %>% melt()
colnames(V4V5_id_taxa)[1] <- c("Metric")
V4V5_id_taxa$Region <- "V4V5"
V4V5_id_taxa$variable <- "IDTAXA"
colnames(V4V5_id_taxa)
V4V5_id_taxa$Dataset <- "Holdout"
V4V5_id_taxa3 <- V4V5_id_taxa %>% filter(!Metric %in% c("Unassigned Correct", "Unassigned Incorrect", "Accuracy"))
V4V5_id_taxa3_full <- V4V5_id_taxa

V4V5_id_taxa <- rbind(V4V5_id_taxa1, V4V5_id_taxa2, V4V5_id_taxa3)

plot_data <- rbind(plot_data, FL_id_taxa, V4V5_id_taxa, V1V2_id_taxa)

cols <- c("Parathaa Specific (default)"="#FA8072",
          "Parathaa Sensitive"="#660000",
          "Naive Bayes-EM"="#143d59",
          "Naive Bayes Multi-EM"="#1e1eff",
          "Naive Bayes"="#47008e",
          "IDTAXA"="#BA8E23")

shapes <- c("Parathaa Specific (default)"=16,
            "Parathaa Sensitive"=16,
            "Naive Bayes-EM"=17,
            "Naive Bayes"=17,
            "Naive Bayes Multi-EM"=17,
            "IDTAXA"=18)


small_plot <- ggplot(plot_data, aes(x=Dataset, y=as.numeric(value), group=variable, color=variable, shape=variable)) + 
  geom_jitter(position=position_jitterdodge(dodge.width = 1, jitter.width = .1), size=3, alpha=0.6) +
  #geom_point(stat="summary", fun="mean", shape="_", position = position_dodge(width = 1), size=10) +
  facet_grid(rows=vars(Metric), cols=vars(Region), scales="free_y") +
  theme_bw() + 
  scale_color_manual(values = cols, name="Software", breaks=
                       c("Parathaa Specific (default)", "Parathaa Sensitive", "Naive Bayes-EM", "Naive Bayes Multi-EM", 
                        "Naive Bayes", "IDTAXA"))+ 
  ylab("") +
  ggtitle("Species") + theme(strip.text.y=element_text(angle=0), axis.text.x = element_text(angle=65, vjust = 0.5)) +
  scale_shape_manual(values=shapes, name="Software", breaks=
                       c("Parathaa Specific (default)", "Parathaa Sensitive", "DADA2", "DADA2 Multi", 
                         "DADA2 Exact Match", "DADA2 Naive Bayes", "IDTAXA"))+
  scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
  theme(panel.spacing.y = unit(.5, "lines")) +
  xlab("") +
  theme(axis.text.x=element_blank(), axis.ticks.x =  element_blank())

ggsave("~/Dropbox_Harvard/hutlab/Jacob/Paper_Submissions/Parathaa/NAR/Figures/Figure2/Figure2_revised_panelA.pdf",
       width=12, height=6)


### supplemental plot
FL_id_taxa_full <- rbind(FL_id_taxa1_full, FL_id_taxa2_full, FL_id_taxa3_full)
V1V2_id_taxa_full <- rbind(V1V2_id_taxa1_full, V1V2_id_taxa2_full, V1V2_id_taxa3_full)
V4V5_id_taxa_full <- rbind(V4V5_id_taxa1_full, V4V5_id_taxa2_full, V4V5_id_taxa3_full)

plot_dada2 <- rbind(species_plots[[1]]$data, FL_id_taxa_full, V1V2_id_taxa_full, V4V5_id_taxa_full)



full_plot <- ggplot(plot_dada2, aes(x=Dataset, y=as.numeric(value), group=variable, color=variable, shape=variable)) + 
  #geom_jitter(position=position_jitterdodge(dodge.width = 1, jitter.width = .1), size=3, alpha=0.7) +
  geom_point(stat="summary", fun="mean", position = position_dodge(width = 1), size=2) +
  facet_grid(rows=vars(Metric), cols=vars(Region), scales="free_y") +
  theme_bw() + 
  scale_color_manual(values = cols, name="Software", breaks=
                       c("Parathaa Specific (default)", "Parathaa Sensitive", "Naive Bayes-EM", "Naive Bayes Multi-EM", 
                         "Naive Bayes", "IDTAXA"))+ 
  ylab("") +
  ggtitle("Species") + theme(strip.text.y=element_text(angle=0), axis.text.x = element_text(angle=65, vjust = 0.5)) +
  scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
  scale_shape_manual(name="Software", values=shapes, breaks=
                       c("Parathaa Specific (default)", "Parathaa Sensitive", "DADA2", "DADA2 Multi", 
                         "DADA2 Exact Match", "DADA2 Naive Bayes")) +
  geom_rect(xmin=0,
            xmax=1.5,
            ymin=0, ymax=1, color=NA, fill="grey", alpha=0.001) +
  geom_rect(xmin=1.5,
            xmax=2.5,
            ymin=0, ymax=1, color=NA, fill="black", alpha=0.01) +
  geom_rect(xmin=2.5,
            xmax=3.5,
            ymin=0, ymax=1, color=NA, fill="grey", alpha=0.001) +
  geom_rect(xmin=3.5,
            xmax=4.5,
            ymin=0, ymax=1, color=NA, fill="black", alpha=0.01) +
  geom_rect(xmin=4.5,
            xmax=6,
            ymin=0, ymax=1, color=NA, fill="grey", alpha=0.001)
full_plot

ggsave("~/Dropbox_Harvard/hutlab/Jacob/Paper_Submissions/Parathaa/NAR/Figures/Figure2/supplemental_bench.pdf",
       width=14, height=8)
