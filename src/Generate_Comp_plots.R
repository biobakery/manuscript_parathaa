require(docopt)

'Usage:
  Generate_Comp_plots.R [--query <query reads> --abund_tab <abundance table in tsv format> --parathaa_assign <parathaa assignments> --dada_db <dada2 formatted genus DB> --dada_db_sp <dada2 formatted species DB> --output <directory for output> --skip <taxa_level_to_skip> --subset <subset data> --allowMult <allow dada2 multi assignments>]
  
Options:
  --query query reads
  --abund_tab abundance table
  --parathaa_assign parathaa assignment file
  --dada_db dada2 formatted genus level database
  --dada_db_sp dada2 formatted species database
  --output output directory
  --subset should we subset the dataset [default: FALSE]
  --allowMulti should dada2 assign taxonomy with multi [default: TRUE]
  --skip taxonomic level to skip if needed [default: ""]
  ]' -> doc

opts <- docopt(doc)

## Load in libraries that are required
library(dada2)
library(seqinr)
library(vegan)
library(dplyr)
library(stringr)
library(microbiomeutilities)
library(phyloseq)
library(cowplot)

color_hex_codes <- c("#cd5645", "#56a959", "#7967d5", "#7bb136", "#c05db6", 
                              "#d39d29", "#7c87c5", "#b9853c", "#4ba98c", "#c75980", 
                              "#96894d")
query_reads <- "~/Repos/Hills_Project/Parathaa_Project/Parathaa_on_real_data/mine_data/v45/dada2_output/rep-seqs/dna-sequences.fasta"
abundance_file <- "~/Repos/Hills_Project/Parathaa_Project/Parathaa_on_real_data/mine_data/v45/dada2_output/abundance_table/feature-table.tsv"
parathaa_file <- "~/Repos/Hills_Project/Parathaa_Project/Real_Datasets_July_22/MINE/v45/parathaa_assignments/taxonomic_assignments.tsv"

query_reads <- opts$query
abundance_file <- opts$abund_tab
parathaa_file <- opts$parathaa_assign
dada_db <- "~/Repos/Hills_Project/Parathaa_Project/Exploration/ASD_Genus_investigation/input/20231215.silva.seed_v138_1.ng.dada.fasta"
dada_db_sp <- "~/Repos/Hills_Project/Parathaa_Project/Exploration/ASD_Genus_investigation/input/20231215_silva.seed_v138_1.ng.dada.sp.fasta"

dada_db <- opts$dada_db
dada_db_sp <- opts$dada_db_sp

skip <- opts$skip
outputDir <- opts$output

#takes in query reads and assigns taxonomy to them using dada2
assign_dada_taxonomy <- function(query_reads, dada_db, dada_db_sp, allow_multi){
  
  
  names <- names(read.fasta(query_reads))
  #assign with naive bayes
  dada2_assignment <- assignTaxonomy(query_reads, dada_db, multithread = T)
  
  #add speices
  dada2_assignment <- addSpecies(dada2_assignment, dada_db_sp, allowMultiple = allow_multi)
  rownames(dada2_assignment) <- names

  tax_dada3 <- as.data.frame(dada2_assignment) %>% 
    rowwise() %>% 
    mutate(Species = if_else(!is.na(Species), 
                             paste( paste(Genus), str_split(Species, "/",simplify = T), collapse =";"), 
                             NA)
    ) %>% 
    as.matrix()
  
  rownames(tax_dada3) <- rownames(dada2_assignment)
  
  return(tax_dada3)
}


generate_phyloseq_object <- function(parathaa_assignments, dada2_assignment, abundance_table){
  
  ## Make DADA2 phyloseq object ##
  TAX_dada <- tax_table(as.matrix(dada2_assignment))
  rownames(TAX_dada) <- paste0(rownames(TAX_dada), "a")
  otumat <- as.matrix(abundance_table) 
  samp_dada <- data.frame(colnames(otumat), rep("DADA2", length(colnames(otumat))))  
  
  colnames(otumat) <- paste0(colnames(otumat), "a")
  rownames(otumat) <- paste0(rownames(otumat), "a")
  OTU_dada <- otu_table(otumat, taxa_are_rows = TRUE)
  
  colnames(samp_dada) <- c("sampleID", "Taxonomy_type")
  rownames(samp_dada) <-  paste0(samp_dada$sampleID, "a")
  SAMP_dada <- sample_data(samp_dada)
  
  ps1_dada <- phyloseq(OTU_dada, TAX_dada, SAMP_dada)
  
  ## Make parathaa phyloseq object ##
  tax_parathaa <- parathaa_assignments %>% select(Kingdom, Phylum, Class, Order, Family, Genus, Species)
  TAX_parathaa <- tax_table(as.matrix(tax_parathaa))
  rownames(TAX_parathaa) <- paste0(rownames(TAX_parathaa), "b")
  otumat <- as.matrix(abundance_table)
  samp_parathaa <- data.frame(colnames(otumat), rep("Parathaa", length(colnames(otumat))))  


  colnames(otumat) <- paste0(colnames(otumat), "b")
  rownames(otumat) <- paste0(rownames(otumat), "b")
  OTU_parathaa <- otu_table(otumat, taxa_are_rows = TRUE)
  
  colnames(samp_parathaa) <- c("sampleID", "Taxonomy_type")
  rownames(samp_parathaa) <- paste0(samp_parathaa$sampleID, "b")  
  SAMP_parathaa <- sample_data(samp_parathaa)
  
  ps1_parathaa <- phyloseq(OTU_parathaa, TAX_parathaa, SAMP_parathaa)  
  
  ps1_all<- merge_phyloseq(ps1_dada, ps1_parathaa)
  
  return(ps1_all)

}

### Generate bray PCoAs at each level
generate_bray_PCoA <- function(ps1.com, skip="", remove_unknown=FALSE){
  
  Levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  plotList <- list()
  for(level in Levels){
    message(level)
    if(skip==level){
      break
    }
    if(remove_unknown){
      ps1.com <- subset_taxa(ps1.com,  eval(as.name(level)) != "unassigned")
    }
    ps1.com.lev <- aggregate_rare(ps1.com, level, detection = .001/100, prevalence = 1/100)
    ps1.com.lev.agg <- aggregate_taxa(ps1.com.lev, level=level)
    ps1.com.lev.agg.rel <- microbiome::transform(ps1.com.lev.agg, "compositional")
    GP.ord <- ordinate(ps1.com.lev.agg.rel, "PCoA", "bray")
    #brays <- phyloseq::distance(ps1.com.lev.agg.rel, method="bray")
    #adonis2(brays ~ sample_data(ps1.com.lev.agg.rel)$Taxonomy_type)
    plot_data = plot_ordination(ps1.com.lev.agg.rel, GP.ord, type="samples", color="Taxonomy_type", title=level, justDF = TRUE) 
    
    p1 <- ggplot(plot_data, aes(x=Axis.1, y=Axis.2, color=Taxonomy_type)) + geom_point(alpha=0.3, size=2) +
      theme_bw() + geom_polygon(aes(group=sampleID), color="grey", alpha=0.1)
    
    cum_eg <- sum(GP.ord$values$Eigenvalues)
    var1 <- round(GP.ord$values$Eigenvalues[1]/cum_eg * 100, digits = 3)
    var2 <- round(GP.ord$values$Eigenvalues[2]/cum_eg *100, digits = 3)
    
    p1 <- p1 + xlab(paste0("PC1 ", var1)) + ylab(paste0("PC2 ", var2)) + ggtitle(level)
    plotList[[level]] <- p1
    
  }
  
  return(plotList)
}

generate_taxonomy_plots <- function(ps1.com, order){
  plotList <- list()
  Levels <- c("Genus", "Species")
  for(level in Levels){
    
    ps1.com.rel <- microbiome::transform(ps1.com, "compositional")
    ps1.com.rel.lev <- aggregate_rare(ps1.com.rel, level, detection = .001/100, prevalence = 1/100)
    ps1.com.rel.lev.agg <- aggregate_top_taxa2(ps1.com.rel.lev, 11, level)
    
    plot.composition.relAbun <- plot_composition(ps1.com.rel.lev.agg,
                                                 otu.sort = "abundance",
                                                 x.label = "Taxonomy_type") 
    
    
    plot.composition.relAbun$data$Sample <- gsub("a$", "", plot.composition.relAbun$data$Sample)
    plot.composition.relAbun$data$Sample <- gsub("b$", "", plot.composition.relAbun$data$Sample)
    
    plot.composition.relAbun$data$Sample <- factor(plot.composition.relAbun$data$Sample, levels=Genus_order$sampleID)
    
    plot.composition.relAbun <- ggplot(plot.composition.relAbun$data, aes(x=Sample, y=Abundance, fill=Tax)) + 
      geom_bar(stat="identity")
    
    xlabs <- rep("", nrow(sample_data(ps1.com.rel.lev.agg)))
    plot.composition.relAbun <- plot.composition.relAbun + theme(legend.position = "bottom") 
    plot.composition.relAbun <- plot.composition.relAbun +  theme_bw()
    plot.composition.relAbun <- plot.composition.relAbun + ggtitle(NULL) + ylab("Relative Abundance") + labs(fill=level) +
      guides(fill= 
               guide_legend(
                 ncol=1,
                 label.theme = element_text(face="italic", size=9), 
                 title.theme = element_text(size=14),
               )
      )
    plot.composition.relAbun$data$Tax <- as.character(plot.composition.relAbun$data$Tax)
    plot.composition.relAbun$data$Tax <- gsub("(;[^;]*);", "\\1\n", plot.composition.relAbun$data$Tax)
    plot.composition.relAbun$data$Tax <- factor(plot.composition.relAbun$data$Tax)
    plot.composition.relAbun$data$Tax <- relevel(plot.composition.relAbun$data$Tax, "Unknown")


    plot.composition.relAbun <- plot.composition.relAbun +   theme(axis.text.x = element_blank()) + 
      scale_x_discrete(labels=xlabs) + xlab("") + facet_grid(rows=vars(xlabel), scales = "free", space="free") + 
      scale_fill_manual(values= c("white", color_hex_codes))
    plotList[[level]] <- plot.composition.relAbun
    
  }
  return(plotList)
  
  
}

plot_average_taxa <- function(ps1.com){
  
  plotList <- list()
  Levels <- c("Genus", "Species")
  for(level in Levels){
    
    ps1.com.rel <- microbiome::transform(ps1.com, "compositional")
    ps1.com.rel.lev <- aggregate_rare(ps1.com.rel, level, detection = .001/100, prevalence = 1/100)
    ps1.com.rel.lev.agg <- aggregate_top_taxa2(ps1.com.rel.lev, 11, level)
    
    plot.composition.relAbun <- plot_composition(ps1.com.rel.lev.agg,
                                                 otu.sort = "abundance",
                                                 x.label = "Taxonomy_type",
                                                 average_by = "Taxonomy_type") 
    
    
    plot.composition.relAbun$data$Sample <- gsub("a$", "", plot.composition.relAbun$data$Sample)
    plot.composition.relAbun$data$Sample <- gsub("b$", "", plot.composition.relAbun$data$Sample)
    
    plot.composition.relAbun$data$Sample <- factor(plot.composition.relAbun$data$Sample, levels=Genus_order$sampleID)
    
    plot.composition.relAbun <- ggplot(plot.composition.relAbun$data, aes(x=xlabel, y=Abundance, fill=Tax)) + 
      geom_bar(stat="identity")
    
    plot.composition.relAbun <- plot.composition.relAbun + theme(legend.position = "bottom") 
    plot.composition.relAbun <- plot.composition.relAbun +  theme_bw()
    plot.composition.relAbun <- plot.composition.relAbun + ggtitle(NULL) + ylab("Relative Abundance") + labs(fill=level) +
      guides(fill= 
               guide_legend(
                 ncol=1,
                 label.theme = element_text(face="italic", size=9), 
                 title.theme = element_text(size=14),
               )
      )
    plot.composition.relAbun$data$Tax <- as.character(plot.composition.relAbun$data$Tax)
    plot.composition.relAbun$data$Tax <- gsub("(;[^;]*);", "\\1\n", plot.composition.relAbun$data$Tax)
    plot.composition.relAbun$data$Tax <- factor(plot.composition.relAbun$data$Tax)
    plot.composition.relAbun$data$Tax <- relevel(plot.composition.relAbun$data$Tax, "Unknown")
    
    
    plot.composition.relAbun <- plot.composition.relAbun + xlab("") +
      scale_fill_manual(values= c("white", color_hex_codes))
    plotList[[level]] <- plot.composition.relAbun
    
  }
  return(plotList)
  
}

get_variance_by_tool <- function(ps1.com, level, remove_unknown=FALSE){
  
  if(remove_unknown){
    ps1.com <- subset_taxa(ps1.com, eval(as.name(level)) != "unassigned")
  }
  metadata <- as(sample_data(ps1.com), "data.frame")
  
  ps1_Genus <- aggregate_taxa(ps1.com, level = level)
  
  res <- adonis2(distance(ps1_Genus, method="bray") ~ Taxonomy_type, data=metadata)
  
}

#function that takes in a string and checks if all levels are unclassified
check_level <- function(x){
  
  #split the string by ;
  x_split <- str_split(x, ";")
  #number of assignments
  num_assignments <- length(x_split[[1]])
  #get number of unclassified assignments
  num_unclass <- length(which(grepl("Unclassified", x_split[[1]])))

  #if all are unclassified then return NA
  if(num_assignments == num_unclass){
    return(NA)
  #else return original label
  }else{
    return(x)
  }
}

#read in Parathaa assignments
parathaa_assignments <- read.table(parathaa_file, sep="\t", header=T, row.names=1, stringsAsFactors = FALSE)
#remove unclassified labels when all things at that level are unclassified
taxa_levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
parathaa_assignments_fix <- parathaa_assignments %>% rowwise %>% mutate_at(vars(taxa_levels), ~ check_level(.))
parathaa_assignments_fix <- data.frame(parathaa_assignments_fix)
rownames(parathaa_assignments_fix) <- rownames(parathaa_assignments)
parathaa_assignments_fix[parathaa_assignments_fix==""] <- NA
parathaa_assignments_fix[parathaa_assignments_fix=="uncultured"] <- NA 
## set unclassified levels to blank..




if(opts$allowMulti=="true" | opts$allowMulti=="True" | opts$allowMulti=="TRUE" | opts$allowMulti=="T"){
  opts$allowMulti <- TRUE
}else
  opts$allowMulti <- FALSE

#get dada2 assignments
dada2_assignment <- assign_dada_taxonomy(query_reads, dada_db, dada_db_sp, opts$allowMulti)
dada2_assignment[dada2_assignment==""] <- NA
dada2_assignment[dada2_assignment=="uncultured"] <- NA

abundance_table <- read.table(abundance_file, sep="\t", header=T, check.names = F, comment.char="", skip=1, row.names=1)



#if we want to subset to 100 samples (used in the oral dataset that is very large)
if(opts$subset==TRUE){
  abundance_table <- abundance_table[,sample(colnames(abundance_table, 100))] 
}


merged_phyloseq <- generate_phyloseq_object(parathaa_assignments_fix, dada2_assignment, abundance_table)


print(merged_phyloseq)


#generate output directory if it doesn't already exist
dir.create(outputDir, recursive = T, showWarnings = F)

#save phyloseq object
saveRDS(merged_phyloseq, file=file.path(outputDir, "merged_phyloseq.RDS"))

bray_plots<- generate_bray_PCoA(merged_phyloseq, skip=skip, remove_unknown = FALSE)

saveRDS(bray_plots, file=file.path(outputDir, "bray_plots.RDS"))

## add sorting info to merge so make barplots look better... (we can sort it by genus abundance)
if(skip=="Species"){
  Genus_order <- bray_plots[[5]]$data %>% filter(grepl("Parathaa", Taxonomy_type)) %>% arrange(Axis.1)
}else{
  Genus_order <- bray_plots[[6]]$data %>% filter(grepl("Parathaa", Taxonomy_type)) %>% arrange(Axis.1)
}

taxa_plots <- generate_taxonomy_plots(merged_phyloseq, Genus_order)
saveRDS(taxa_plots, file=file.path(outputDir, "taxa_bars.RDS"))

average_plots <- plot_average_taxa(merged_phyloseq)
saveRDS(average_plots, file=file.path(outputDir, "average_plots.RDS"))

## done saved all the plots

if(opts$subset){
  set.seed(29)
  abundance_table <- abundance_table[,sample(1:ncol(abundance_table), 100)]
  
  merged_phyloseq <- generate_phyloseq_object(parathaa_assignments_fix, dada2_assignment, abundance_table)
  
  
  print(merged_phyloseq)
  
  #save phyloseq object
  saveRDS(merged_phyloseq, file=file.path(outputDir, "merged_phyloseq_sub.RDS"))
  
  bray_plots<- generate_bray_PCoA(merged_phyloseq, skip=skip, remove_unknown = FALSE)
  
  saveRDS(bray_plots, file=file.path(outputDir, "bray_plots_sub.RDS"))
  
  ## add sorting info to merge so make barplots look better... (we can sort it by genus abundance)
  if(skip=="Species"){
    Genus_order <- bray_plots[[5]]$data %>% filter(grepl("Parathaa", Taxonomy_type)) %>% arrange(Axis.1)
  }else{
    Genus_order <- bray_plots[[6]]$data %>% filter(grepl("Parathaa", Taxonomy_type)) %>% arrange(Axis.1)
  }
  
  taxa_plots <- generate_taxonomy_plots(merged_phyloseq, Genus_order)
  saveRDS(taxa_plots, file=file.path(outputDir, "taxa_bars_sub.RDS"))
  
  average_plots <- plot_average_taxa(merged_phyloseq)
  saveRDS(average_plots, file=file.path(outputDir, "average_plots_sub.RDS"))
  
}


# ps1_dada_Genus <- aggregate_taxa(ps1_dada, level="Genus")
# ps1_dada_Genus <- transform(ps1_dada_Genus, transform = "compositional")
# ps1_genus_rowmeans <- rowMeans(ps1_dada_Genus@otu_table)
# 
# 
# ps1_parathaa_Genus <- aggregate_taxa(ps1_parathaa, level="Genus")
# ps1_parathaa_Genus <- transform(ps1_parathaa_Genus, transform = "compositional")
# ps1_parathaa_rowmeans <- rowMeans(ps1_parathaa_Genus@otu_table)
