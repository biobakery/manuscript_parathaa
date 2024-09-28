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

#opts$specific <- "~/Repos/Hills_Project/Parathaa_Project/Main_results_July_22/output_specific/benchmarks/"
#opts$sensitive <- "~/Repos/Hills_Project/Parathaa_Project/Main_results_July_22/output/benchmarks/"

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

## function that plots genus benchmarks
plot_genus_figures <- function(genus_bench_list){
  
  #load in genus benchmark data
  full_df_melt <- melt_bench_data(genus_bench_list)
  
  #rename the datasets to be more legible
  full_df_melt$Dataset <- gsub("holdoutOG", "Orig. Holdout", full_df_melt$Dataset)
  full_df_melt$Dataset <- gsub("holdout[1-3]", "Holdout", full_df_melt$Dataset)
  full_df_melt$Dataset <- gsub("original", "Original", full_df_melt$Dataset)
  full_df_melt$Dataset <- gsub("novel", "Genera outside of DB", full_df_melt$Dataset)
  full_df_melt$Dataset <- gsub("even", "Single represenative genera", full_df_melt$Dataset)
  
  #filter out one of the Parathaa measurements as its repeat information
  full_df_melt <- full_df_melt %>% filter(Region != "Full_naive" | variable != "Parathaa")
  
  
  full_df_melt$variable <- as.character(full_df_melt$variable)
  #rename the DADA2 benchmarks to their true names for Genus level
  full_df_melt$variable[which(full_df_melt$Region=="Full_naive" & full_df_melt$variable=="DADA2")] <- "DADA2 minboot=80"
  full_df_melt$variable[which(full_df_melt$Region=="Full_exact" & full_df_melt$variable=="DADA2")] <- "DADA2 minboot=50"
  
  #clean region names
  full_df_melt$Region <- gsub("_.*", "", full_df_melt$Region)
  

  full_df_melt$value <- as.numeric(full_df_melt$value)
  full_df_melt$Metric <- factor(full_df_melt$Metric, 
                                levels = c("F1 Score", "Precision", "Recall", "Accuracy", "Uniquely Correct", 
                                           "One-to-many Correct", "Incorrect", "Unassigned Correct", "Unassigned Incorrect"))
  
  full_df_melt$Dataset <- factor(full_df_melt$Dataset, levels=c("Original", "Orig. Holdout", "Genera outside of DB",
                                                                "Single represenative genera", "Holdout"))
  
  
  full_plot <- ggplot(full_df_melt, aes(x=Dataset, y=as.numeric(value), group=variable, color=variable, shape=variable)) + 
    geom_jitter(position=position_jitterdodge(dodge.width = 1, jitter.width = .1), size=3, alpha=1) +
    #geom_point(stat="summary", fun="mean", shape="_", position = position_dodge(width = 1), size=10) +
    facet_grid(rows=vars(Metric), cols=vars(Region), scales="free_y") +
    theme_bw() + 
    scale_color_manual(name="Software", values=c(Parathaa="#A41034", DADA2="#143d59", "DADA2 minboot=50"="#3F054F", "DADA2 minboot=80"="#2a7fb9"))+ 
    ylab("") +
    ggtitle("Genus") + theme(strip.text.y=element_text(angle=0), axis.text.x = element_text(angle=65, vjust = 0.5)) +
    scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
    scale_shape_manual(name="Software", values=c(DADA2=17, "DADA2 minboot=50"=17, "DADA2 minboot=80"=17, "Parathaa"=16)) +
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
  
  

  
  smaller_plot <- ggplot(filt_data, aes(x=Dataset, y=as.numeric(value), group=variable, color=variable, shape=variable)) + 
    geom_jitter(position=position_jitterdodge(dodge.width = 1, jitter.width = .1), size=3, alpha=0.7) +
    #geom_point(stat="summary", fun="mean", shape="_", position = position_dodge(width = 1), size=10) +
    facet_grid(rows=vars(Metric), cols=vars(Region), scales="free_y") +
    theme_bw() + 
    scale_color_manual(name="Software", values=c(Parathaa="#A41034", DADA2="#143d59", "DADA2 minboot=50"="#3F054F", "DADA2 minboot=80"="#2a7fb9"))+ 
    ylab("") +
    ggtitle("Genus") + theme(strip.text.y=element_text(angle=0), axis.text.x = element_text(angle=65, vjust = 0.5)) +
    scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
    scale_shape_manual(name="Software", values=c(DADA2=17, "DADA2 minboot=50"=17, "DADA2 minboot=80"=17, "Parathaa"=16)) +
    geom_blank(data=dumby_data)
  
  return(list(full_plot, smaller_plot))

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
  
  
  
  cols <- c("Parathaa Specific (default)"="#A41034",
            "Parathaa Sensitive"="#FFA500",
            "DADA2"="#143d59",
            "DADA2 Multi"="#3F054F",
            "DADA2 Naive Bayes"="#ca5cdd",
            "DADA2 Exact Match"="#82EEFD")
  
  shapes <- c("Parathaa Specific (default)"=16,
              "Parathaa Sensitive"=16,
              "DADA2"=17,
              "DADA2 Multi"=17,
              "DADA2 Naive Bayes"=17,
              "DADA2 Exact Match"=17)
  
  full_df_melt$Dataset <- gsub("holdout[1-3]", "Holdout", full_df_melt$Dataset)
  full_df_melt$Dataset <- gsub("original", "Original", full_df_melt$Dataset)
  full_df_melt$Dataset <- gsub("holdoutOG", "Orig. Holdout", full_df_melt$Dataset)
  full_df_melt$Dataset <- gsub("novel", "Genera outside of DB", full_df_melt$Dataset)
  full_df_melt$Dataset <- gsub("even", "Single represenative genera", full_df_melt$Dataset)
  
  full_df_melt$variable <- factor(full_df_melt$variable, levels=c("Parathaa Specific (default)", 
                                                                  "DADA2 Multi",
                                                                  "DADA2 Exact Match",
                                                                  "Parathaa Sensitive",
                                                                  "DADA2",
                                                                  "DADA2 Naive Bayes"
                                                                  ))
  
  full_df_melt$Dataset <- factor(full_df_melt$Dataset, levels=c("Original", "Orig. Holdout", "Holdout",
                                                                "Single represenative genera",  "Genera outside of DB"))
  
  full_plot <- ggplot(full_df_melt, aes(x=Dataset, y=as.numeric(value), group=variable, color=variable, shape=variable)) + 
    geom_jitter(position=position_jitterdodge(dodge.width = 1, jitter.width = .1), size=3, alpha=0.7) +
    #geom_point(stat="summary", fun="mean", shape="_", position = position_dodge(width = 1), size=10) +
    facet_grid(rows=vars(Metric), cols=vars(Region), scales="free_y") +
    theme_bw() + 
    scale_color_manual(values = cols, name="Software")+ 
    ylab("") +
    ggtitle("Species") + theme(strip.text.y=element_text(angle=0), axis.text.x = element_text(angle=65, vjust = 0.5)) +
    scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
    scale_shape_manual(name="Software", values=shapes) +
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
    Metric %in% top_panel & variable == "Parathaa" ~ 1,
    Metric %in% top_panel & variable == "DADA2" ~ 0.5,
    Metric %in% top_panel & variable == "DADA2 minboot=50" ~ 0.5,
    Metric %in% bottom_panel & variable == "Parathaa" ~ 0.5,
    Metric %in% bottom_panel & variable == "DADA2" ~ 0.5,
    Metric %in% bottom_panel & variable == "DADA2 minboot=50" ~ 0.5,
  ))
  
  small_plot <- ggplot(filt_data, aes(x=Dataset, y=as.numeric(value), group=variable, color=variable, shape=variable)) + 
    geom_jitter(position=position_jitterdodge(dodge.width = 1, jitter.width = .2), size=3, alpha=0.6) +
    #geom_point(stat="summary", fun="mean", shape="_", position = position_dodge(width = 1), size=10) +
    facet_grid(rows=vars(Metric), cols=vars(Region), scales="free_y") +
    theme_bw() + 
    scale_color_manual(values = cols, name="Software")+ 
    ylab("") +
    ggtitle("Species") + theme(strip.text.y=element_text(angle=0), axis.text.x = element_text(angle=65, vjust = 0.5)) +
    scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
    scale_shape_manual(values=shapes, name="Software")+
    geom_blank(data=dumby_data)
  
  
  return(list(full_plot, small_plot))
  
}


mult_adjusted_Genus <- get_data(opts$specific, level="Genus", adjust_name = "adjust", mult = "mult")
mult_adjusted_Species <- get_data(opts$specific, level="Species", adjust_name="adjust", mult="mult")
nomult_adjusted_Species <- get_data(opts$specific, level="Species", adjust_name = "adjust", mult="nomult")
sensitive_Species <- get_data(opts$sensitive, level="Species", adjust_name = "adjust", mult="mult")


genus_plots <- plot_genus_figures(mult_adjusted_Genus)

##save genus plots
ggsave(filename = file.path(opts$output, "Genus_full_data.pdf"), plot=genus_plots[[1]], height = 8.25, width=1, units = "in")
ggsave(filename = file.path(opts$output, "Genus_filt_data.pdf"), plot=genus_plots[[2]], height = 6, width=8, units = "in")


species_plots <- plot_benchs_species(mult_bench = mult_adjusted_Species,
                                     no_mult_bench = nomult_adjusted_Species,
                                     sens_bench = sensitive_Species)

ggsave(filename = file.path(opts$output, "Species_full_data.pdf"), plot=species_plots[[1]], height = 8.25, width=11, units = "in")
ggsave(filename = file.path(opts$output, "Species_filt_data.pdf"), plot=species_plots[[2]], height = 6, width=8, units = "in")


