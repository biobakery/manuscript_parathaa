## Compare taxonomies from SILVA vs PARATHAA

require(docopt)

'Usage:
  Plots.Figure.3.R [-o <output> --counts <closed_ref_taxonomy_table> --fasta <fasta_file> --dada_db <dada2 db> --dada_db_sp <dada2 spec db> --paraAssign <parathaa_assignments> --countSeq <count file with seqs>]
  
  
Options:
  --counts load all_samples_taxonomy_closed_reference.tsv 
  --countSeq count file with sequences
  --fasta load amplicons_fromKelsey.fasta (fasta file for the counts table)
  --dada_db dada2 db
  --dada_db_sp dada2 Species db
  --paraAssign parathaa assignments
  -o output
  ]' -> doc


opts <- docopt(doc)


library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(dplyr) # data handling  
library(treeio)
library(gridExtra)
library(lemon)
library(vegan)
library(microbiomeutilities)
library(pheatmap)
library(dada2)
library(devtools)
library(pals)
library(ComplexHeatmap)
library(stringr)
library(tidyr)
suppressPackageStartupMessages(library(seqinr))
source_url("https://raw.githubusercontent.com/lrjoshi/FastaTabular/master/fasta_and_tabular.R")

#### SET OUTPUT DIRECTORY FOR PLOTS ####
outputDir <- opts$o#"output/Figures/ASD/allowMult"
dadaAllowMult <- TRUE

############################
### ASD Dataset Analysis ###
############################

## Read in data with taxonomy from Kelsey's analysis with SILVA
kelseyDataRaw <- read.delim(opts$counts, #"./input/all_samples_taxonomy_closed_reference.tsv", 
                            sep='\t', fill=T, stringsAsFactors = F, header=T)

message("read in counts")

tax_oldDADA2 <- kelseyDataRaw %>%
  dplyr::select(taxonomy) %>%
  separate(taxonomy, into=c("Kingdom", "Phylum", "Class", "Order", "Family",
                            "Genus", "Species"), sep=";")

## Remove samples as instructed by Kelsey
kelseyDataFull <- kelseyDataRaw %>%
  dplyr::select(-c(X,
                   CBSF036_S66_L001, 
                   CBSF060_S82_L001, 
                   CBSF078_S77_L001, 
                   CBSF079_S65_L001,
                   Undetermined_S0_L001,
                   taxonomy))

## SILVA Seed-based DADA2 training set
taxa <- assignTaxonomy(opts$fasta,#"./input/amplicons_fromKelsey.fasta", 
                       opts$dada_db, #"./input/20231215.silva.seed_v138_1.ng.dada.fasta",
                       multithread=TRUE)

taxa.sp <- addSpecies(taxa, opts$dada_db_sp, #"./input/20231215_silva.seed_v138_1.ng.dada.sp.fasta", 
                      allowMultiple = dadaAllowMult)

tax_dada <- as.matrix(taxa.sp)

## Get taxa IDs
kelseyData <- read.delim(opts$countSeq, #"input/all_samples_taxonomy_closed_reference_withseqs.tsv", 
                         sep='\t', fill=T, stringsAsFactors = F)

message("read in count Seq")

kelseyData <- kelseyData[c(1,126)] 

taxaIDs <- paste("ID", 1:nrow(kelseyData), sep = "")
kelseyData <- cbind(taxaIDs, kelseyData[,2], kelseyData[,1])
kelseyData <- as.data.frame(kelseyData)
colnames(kelseyData) <- c("taxaIDs", "Taxonomy_1", "seq")

#Get IDs from format_amplicon_table.R "kelseyData" object
rownames(tax_dada) <- paste0(kelseyData$taxaIDs, "a")

tax_dada <- cbind(tax_dada, "Species2"=NA)
tax_dada[which(!is.na(tax_dada[,"Species"])), "Species2"] <-  
  paste(tax_dada[which(!is.na(tax_dada[,"Species"])), "Genus"], 
        tax_dada[which(!is.na(tax_dada[,"Species"])), "Species"] )
tax_dada <- tax_dada[,!colnames(tax_dada) %in% "Species"]
colnames(tax_dada)[which(colnames(tax_dada)=="Species2")] <- "Species"

TAX_dada <- tax_table(tax_dada)


otumat <- as.matrix(kelseyDataFull)
rownames(otumat) <- paste0(kelseyData$taxaIDs, "a")

samp_dada <- data.frame(colnames(otumat), rep("DADA2", length(colnames(otumat))))

colnames(otumat) <- paste0(colnames(otumat), "a")
OTU_dada <- otu_table(otumat, taxa_are_rows = TRUE)


colnames(samp_dada) <- c("sampleID", "Taxonomy_type")
rownames(samp_dada) <-  paste0(samp_dada$sampleID, "a")
str(samp_dada)
SAMP_dada <- sample_data(samp_dada)

ps1_dada <- phyloseq(OTU_dada, TAX_dada, SAMP_dada)

print(ps1_dada)
#save(ps1_dada, file="output/ps1_dada.RData")
#load("output/ps1_dada.RData")

## Read in PARATHAA results
parathaData <- read.delim(opts$paraAssign, #"./output/20231203_kindom_fix_ASD/taxonomic_assignments.tsv", 
                          sep='\t', fill=T, stringsAsFactors = F, header=T)

message("read in parathaa results")

parathaData$query.num <- as.numeric(gsub("ID", "", parathaData$query.name))
parathaData <- parathaData %>% arrange(query.num)

tax_paratha <- parathaData %>%
  dplyr::select(query.name, Kingdom, Phylum, Class, Order, Family, Genus, Species) 

taxmat_paratha <- as.matrix(tax_paratha)
rownames(taxmat_paratha) <- paste0(tax_paratha$query.name, "b")
taxmat_paratha <- taxmat_paratha[,-1]
TAX_parathaa <- tax_table(taxmat_paratha)

otumat_parathaa <- kelseyDataFull

otumat_parathaa <- as.matrix(otumat_parathaa)
samp_parathaa <- data.frame(colnames(otumat_parathaa), rep("parathaa", length(colnames(otumat_parathaa))))
rownames(otumat_parathaa) <- paste0(kelseyData$taxaIDs, "b")
colnames(otumat_parathaa) <- paste0(colnames(otumat_parathaa), "b")
OTU_parathaa <- otu_table(otumat_parathaa, taxa_are_rows = TRUE)

colnames(samp_parathaa) <- c("sampleID", "Taxonomy_type")
rownames(samp_parathaa) <- paste0(samp_parathaa$sampleID, "b")
str(samp_parathaa)
SAMP_parathaa <- sample_data(samp_parathaa)

ps1_parathaa <- phyloseq(OTU_parathaa, TAX_parathaa, SAMP_parathaa)


ps1_all<- merge_phyloseq(ps1_dada, ps1_parathaa)
ps1.com <- ps1_all

# We need to set Palette
taxic <- as.data.frame(ps1.com@tax_table) # this will help in setting large color options

taxic$OTU <- rownames(taxic) # Add the OTU ids from OTU table into the taxa table at the end.
colnames(taxic) # You can see that we now have extra taxonomy levels.

taxmat <- as.matrix(taxic) # convert it into a matrix.
new.tax <- tax_table(taxmat) # convert into phyloseq compatible file.
tax_table(ps1.com) <- new.tax # incroporate into phyloseq Object

# now edit the unclassified taxa
tax_table(ps1.com)[is.na(tax_table(ps1.com)[, "Genus"]), "Genus"] <- "Unknown"
tax_table(ps1.com)[is.na(tax_table(ps1.com)[, "Family"]), "Family"] <- "Unknown"
tax_table(ps1.com)[is.na(tax_table(ps1.com)[, "Order"]), "Order"] <- "Unknown"
tax_table(ps1.com)[is.na(tax_table(ps1.com)[, "Class"]), "Class"] <- "Unknown"
tax_table(ps1.com)[is.na(tax_table(ps1.com)[, "Phylum"]), "Phylum"] <- "Unknown"
tax_table(ps1.com)[is.na(tax_table(ps1.com)[, "OTU"]), "OTU"] <- "Unknown"
tax_table(ps1.com)[is.na(tax_table(ps1.com)[, "Species"]), "Species"] <- "Unknown"
tax_table(ps1.com)[is.na(tax_table(ps1.com)[, "Kingdom"]), "Kingdom"] <- "Unknown"


# it would be nice to have the Taxonomic names in italics.
# for that we set this

guide_italics <- guides(fill = guide_legend(label.theme = element_text(
  size = 15,
  face = "italic", colour = "Black", angle = 0
)))


## Species and Genus level plots
# first remove the phy_tree

ps1.com@phy_tree <- NULL

Levels <- c("Genus", "Species")
plotList <- list()
# Taxonomy plots
for(level in Levels){
  
  ps1.com.rel <- microbiome::transform(ps1.com, "compositional")
  ps1.com.rel.lev <- aggregate_rare(ps1.com.rel, level, detection = .001/100, prevalence = 1/100)
  ps1.com.rel.lev.agg <- aggregate_top_taxa2(ps1.com.rel.lev, 11, level)
  
  plot.composition.relAbun <- plot_composition(ps1.com.rel.lev.agg,
                                               sample.sort="Taxonomy_type",
                                               otu.sort = "abundance",
                                               x.label = "sampleID") 
  xlabs <- rep("", nrow(sample_data(ps1.com.rel.lev.agg)))
  xlabs[ceiling(nrow(sample_data(ps1.com.rel.lev.agg))/4)] <- "DADA2"
  xlabs[ceiling(nrow(sample_data(ps1.com.rel.lev.agg))*3/4)+1] <- "parathaa"
  plot.composition.relAbun <- plot.composition.relAbun + theme(legend.position = "bottom") 
  plot.composition.relAbun <- plot.composition.relAbun +  theme_bw()  +  scale_fill_brewer(palette="Paired", labels = ~ str_wrap(.x, width = 50, whitespace_only = F))  
  plot.composition.relAbun <- plot.composition.relAbun + ggtitle(NULL) + ylab("Relative Abundance") + labs(fill=level) +
    guides(fill= 
             guide_legend(
               ncol=1,
               label.theme = element_text(face="italic", size=10), 
               title.theme = element_text(size=14),
             )
    )
  plot.composition.relAbun <- plot.composition.relAbun +   theme(axis.text.x = element_text(size=12, angle=0, hjust=1)) + scale_x_discrete(labels=xlabs) 
  if(level!="Species")
    ggsave(filename = file.path(outputDir, paste0(level, "taxonomy.png")), plot.composition.relAbun, width=8.5, height=6, dpi = 600)
  if(level=="Species")
    ggsave(filename= file.path(outputDir, paste0(level, "taxonomy.png")), plot.composition.relAbun, width=10, height=6, dpi = 600)
  plotList[[level]] <- plot.composition.relAbun
  
}
p1 <- ggarrange(plotList[["Genus"]], plotList[["Species"]], ncol=1, nrow=2, common.legend = F, legend="right", align="hv")

ggsave(filename= file.path(outputDir, "both_taxonomy.png"), p1, width=8, height=8, dpi = 600)

### Ordination plots
Levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
plotList <- list()
for(level in Levels){
  ps1.com.lev <- aggregate_rare(ps1.com, level, detection = .001/100, prevalence = 1/100)
  ps1.com.lev.agg <- aggregate_taxa(ps1.com.lev, level=level)
  ps1.com.lev.agg.rel <- microbiome::transform(ps1.com.lev.agg, "compositional")
  GP.ord <- ordinate(ps1.com.lev.agg.rel, "PCoA", "bray")
  #brays <- phyloseq::distance(ps1.com.lev.agg.rel, method="bray")
  #adonis2(brays ~ sample_data(ps1.com.lev.agg.rel)$Taxonomy_type)
  p1 = plot_ordination(ps1.com.lev.agg.rel, GP.ord, type="samples", color="Taxonomy_type", title=level) + 
    geom_polygon(aes(group=sampleID), color="grey") + theme_bw()
  
  plotList[[level]] <- p1
  
  dad <- otu_table(ps1.com.rel.lev)[,1:119]
  par <- otu_table(ps1.com.rel.lev)[,120:238]
  print("Parathaa unknowns:")
  print(sum(par["Unknown",])/ncol(par))
  print("DADA2 unknowns:")
  print(sum(dad["Unknown",])/ncol(dad))
}

p1 <- ggarrange(plotList[["Phylum"]]+theme(legend.position='hidden'), plotList[["Class"]]+theme(legend.position='hidden'), 
                plotList[["Order"]]+theme(legend.position='hidden'), plotList[["Family"]]+theme(legend.position='hidden'),
                plotList[["Genus"]]+theme(legend.position='hidden'), plotList[["Species"]], ncol=2, nrow=3, common.legend = TRUE, legend="bottom") 

ggsave(filename=file.path(outputDir, "BrayCurtisPCoA.png"), plot=p1, height = 8, width=5, dpi=600)


