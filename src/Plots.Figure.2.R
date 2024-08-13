require(docopt)

'Usage:
  Plots.Figure.2.R [--dada_db <dada2_db> --dada_db_sp <dada2_db_sp> -o <output> --paraAssignV4V5 <V4V5_taxonomy_file> --paraAssignV1V2 <V1V2_taxonomy_file> --fastaV4V5 <V4V5 seqs> --fastaV1V2 <V1V2 seqs> --fastaV1V2uni <unique_V1V2_file> --V1V2Counts <v1v2_count_table> ]
  
  
Options:
  --dada_db location of DADA2 db
  --dada_db_sp location of DADA2 species db
  --fastaV4V5 V4V5 fasta file to assign taxonomy to
  --fastaV1V2uni V1V2 fasta file to assign taxonomy to
  --fastaV1V2 V1V2 fasta file
  --V1V2Counts Counts file for V1V2
  --paraAssignV1V2 parathaa assignments for V1V2
  --paraAssignV4V5 parathaa assignments for V4V5
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


### DADA2 db ###
#DADAdb <- "~/Repos/Hills_Project/Parathaa_Project/Benchmarking/Output_Jun_6/input/20231215.silva.seed_v138_1.ng.dada.fasta"
#DADAdb.sp <- "~/Repos/Hills_Project/Parathaa_Project/Benchmarking/Output_Jun_6/input/20231215_silva.seed_v138_1.ng.dada.sp.fasta"
DADAdb <- opts$dada_db
DADAdb.sp <- opts$dada_db_sp

## Output for plots ##
outputDir <- opts$o

includeSpingo <- FALSE

## utility to refine species names (things that are the same but coded differently by Parathaa/Dada2/SPINGO):
refine_species_names <-function(x){
  x <- x %>% mutate(Species = recode(Species,
                                     "Propionibacterium acnes" = "Cutibacterium acnes",
                                     "Clostridium beijerinckii" = "Clostridium beijerinckii/diolis",
                                     "Clostridium diolis" = "Clostridium beijerinckii/diolis",
                                     "Clostridium sensu stricto 1 diolis" = "Clostridium beijerinckii/diolis",
                                     "Clostridium sensu stricto 1 beijerinckii" = "Clostridium beijerinckii/diolis"
  ))
  return(x)
}

###############################
### Mock community analysis ###
###############################


## Assign taxonomy with dada2 and seed SILVA db
## V4V5
taxaV4V5 <- assignTaxonomy(opts$fastaV4V5, 
                           DADAdb,
                           multithread=TRUE)

nChars <- grep("N", rownames(taxaV4V5))
print(paste("Removing", length(nChars), "sequences with N bases"))
taxaV4V5 <- taxaV4V5[-nChars,]

#Split for memory issues
taxaV4V5.0 <- taxaV4V5[1:21000,]
taxaV4V5.1 <- taxaV4V5[21001:nrow(taxaV4V5),]

taxaV4V5.0.sp <- addSpecies(taxaV4V5.0, DADAdb.sp)
taxaV4V5.1.sp <- addSpecies(taxaV4V5.1, DADAdb.sp)

taxaV4V5.sp <- rbind(taxaV4V5.0.sp, taxaV4V5.1.sp)

taxaV4V5.sp <- cbind(taxaV4V5.sp, "Species2"=NA)
taxaV4V5.sp[which(!is.na(taxaV4V5.sp[,"Species"])), "Species2"] <-  
  paste(taxaV4V5.sp[which(!is.na(taxaV4V5.sp[,"Species"])), "Genus"], 
        taxaV4V5.sp[which(!is.na(taxaV4V5.sp[,"Species"])), "Species"] )
taxaV4V5.sp <- taxaV4V5.sp[,!colnames(taxaV4V5.sp) %in% "Species"]
colnames(taxaV4V5.sp)[which(colnames(taxaV4V5.sp)=="Species2")] <- "Species"

taxaV4V5.sp2 <- refine_species_names(as_tibble(taxaV4V5.sp))
taxaV4V5.sp2 <- as.matrix(taxaV4V5.sp)

otutabV4V5 <- cbind(rownames(taxaV4V5.sp), taxaV4V5.sp)
otutabV4V5 <- as_tibble(otutabV4V5) %>% dplyr::rename("ASV"="V1") %>% dplyr::count(ASV, name="SRR3225703.V4V5") 
otumatV4V5 <- as.matrix(otutabV4V5$SRR3225703.V4V5)
rownames(otumatV4V5) <- otutabV4V5$ASV
colnames(otumatV4V5) <- "SRR3225703.V4V5"

OTU_dada.V4V5 <- otu_table(otumatV4V5, taxa_are_rows = TRUE)

taxtabV4V5 <- cbind(rownames(taxaV4V5.sp), taxaV4V5.sp)
taxtabV4V5 <- as_tibble(taxtabV4V5) %>% 
  dplyr::rename("ASV"="V1") %>% 
  dplyr::group_by(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  dplyr::count(ASV, name="SRR3225703.V4V5") 
taxmatV4V5 <- taxtabV4V5 %>% as.matrix
rownames(taxmatV4V5) <- taxtabV4V5$ASV
taxmatV4V5 <- taxmatV4V5[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]  



TAX_dada.V4V5 <- tax_table(taxmatV4V5)

samp_dada.V4V5 <- data.frame(colnames(OTU_dada.V4V5), 
                             rep("DADA2", length(colnames(OTU_dada.V4V5))),
                             rep("V4V5", length(colnames(OTU_dada.V4V5))))

colnames(samp_dada.V4V5) <- c("sampleID", "Taxonomy_type", "Region")
rownames(samp_dada.V4V5) <- samp_dada.V4V5$sampleID
str(samp_dada.V4V5)
SAMP_dada.V4V5 <- sample_data(samp_dada.V4V5)

ps1_dada.V4V5 <- phyloseq(OTU_dada.V4V5, TAX_dada.V4V5, SAMP_dada.V4V5)

#save(ps1_dada.V4V5, file="output/Figures/mock_analysis/ps1_dada.V4V5.20231205.RData") 

## V1V2 -- in two parts for memory reasons
taxaV1V2.dedup <- assignTaxonomy(opts$fastaV1V2uni, #"./input/SRR3225701.unique.fasta", 
                                 DADAdb,
                                 multithread=TRUE)
nChars <- grep("N", rownames(taxaV1V2.dedup))
print(paste("Removing", length(nChars), "sequences with N bases"))
#taxaV1V2.0 <- taxaV1V2.0[-nChars,]

taxav1v2.dedup.0 <- taxaV1V2.dedup[1:10000,]
taxav1v2.dedup.1 <- taxaV1V2.dedup[10001:20000,]
taxav1v2.dedup.2 <- taxaV1V2.dedup[20001:30000,]
taxav1v2.dedup.3 <- taxaV1V2.dedup[30001:40000,]
taxav1v2.dedup.4 <- taxaV1V2.dedup[40001:50000,]
taxav1v2.dedup.5 <- taxaV1V2.dedup[50001:60000,]
taxav1v2.dedup.6 <- taxaV1V2.dedup[60001:70000,]
taxav1v2.dedup.7 <- taxaV1V2.dedup[70001:80000,]
taxav1v2.dedup.8 <- taxaV1V2.dedup[80001:90000,]
taxav1v2.dedup.9 <- taxaV1V2.dedup[90001:nrow(taxaV1V2.dedup),]


taxaV1V2.0.sp <- addSpecies(taxav1v2.dedup.0, DADAdb.sp)
taxaV1V2.1.sp <- addSpecies(taxav1v2.dedup.1, DADAdb.sp)
taxaV1V2.2.sp <- addSpecies(taxav1v2.dedup.2, DADAdb.sp)
taxaV1V2.3.sp <- addSpecies(taxav1v2.dedup.3, DADAdb.sp)
taxaV1V2.4.sp <- addSpecies(taxav1v2.dedup.4, DADAdb.sp)
taxaV1V2.5.sp <- addSpecies(taxav1v2.dedup.5, DADAdb.sp)
taxaV1V2.6.sp <- addSpecies(taxav1v2.dedup.6, DADAdb.sp)
taxaV1V2.7.sp <- addSpecies(taxav1v2.dedup.7, DADAdb.sp)
taxaV1V2.8.sp <- addSpecies(taxav1v2.dedup.8, DADAdb.sp)
taxaV1V2.9.sp <- addSpecies(taxav1v2.dedup.9, DADAdb.sp)


taxaV1V2.sp <- rbind(taxaV1V2.0.sp, taxaV1V2.1.sp, taxaV1V2.2.sp,
                         taxaV1V2.3.sp, taxaV1V2.4.sp, taxaV1V2.5.sp,
                         taxaV1V2.6.sp,taxaV1V2.7.sp, taxaV1V2.8.sp, taxaV1V2.9.sp)

parsed = seqinr::read.fasta(file(opts$fastaV1V2uni), as.string = TRUE,
                            forceDNAtolower = FALSE, whole.header = FALSE)


table = data.frame("ASV"=unlist(parsed), "query.name" = sapply(parsed, attr, 'name'), row.names=NULL)

v1v2.counts <- read.table(opts$V1V2Counts, #"input/SRR3225701.count_table", 
                          header = T)

otutabV1V2 <- cbind(rownames(taxaV1V2.sp), taxaV1V2.sp)
otutabV1V2 <- otutabV1V2 %>% as.data.frame() %>% dplyr::rename("ASV"="V1") 
otutabV1V2 <- left_join(otutabV1V2, table, by="ASV")
otutabV1V2 <- left_join(otutabV1V2, v1v2.counts, by=c("query.name"="Representative_Sequence")) %>% 
  rename("SRR3225701.V1V2" = "total")

otumatV1V2 <- as.matrix(otutabV1V2$SRR3225701.V1V2)
rownames(otumatV1V2) <- otutabV1V2$ASV
colnames(otumatV1V2) <- "SRR3225701.V1V2"

OTU_dada.V1V2 <- otu_table(otumatV1V2, taxa_are_rows = TRUE)


taxtabV1V2 <- otutabV1V2
taxtabV1V2 <- cbind(taxtabV1V2, "Species2"=NA)
taxtabV1V2[which(!is.na(taxtabV1V2[,"Species"])), "Species2"] <-  
  paste(taxtabV1V2[which(!is.na(taxtabV1V2[,"Species"])), "Genus"], 
        taxtabV1V2[which(!is.na(taxtabV1V2[,"Species"])), "Species"] )
taxtabV1V2 <- taxtabV1V2[,!colnames(taxtabV1V2) %in% "Species"]
colnames(taxtabV1V2)[which(colnames(taxtabV1V2)=="Species2")] <- "Species"
taxtabV1V2 <- refine_species_names(taxtabV1V2)


taxmatV1V2 <- taxtabV1V2 %>% as.matrix
rownames(taxmatV1V2) <- taxtabV1V2$ASV
taxmatV1V2 <- taxmatV1V2[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]  

TAX_dada.V1V2 <- tax_table(taxmatV1V2)

samp_dada.V1V2 <- data.frame(colnames(OTU_dada.V1V2), 
                             rep("DADA2", length(colnames(OTU_dada.V1V2))),
                             rep("V1V2", length(colnames(OTU_dada.V1V2))))
colnames(samp_dada.V1V2) <- c("sampleID", "Taxonomy_type", "Region")
rownames(samp_dada.V1V2) <- samp_dada.V1V2$sampleID
str(samp_dada.V1V2)
SAMP_dada.V1V2 <- sample_data(samp_dada.V1V2)

ps1_dada.V1V2 <- phyloseq(OTU_dada.V1V2, TAX_dada.V1V2, SAMP_dada.V1V2)


## Read in parathaa data

parathaData.V4V5 <- read.delim(opts$paraAssignV4V5, #"./output/20231205_MockV4V5/taxonomic_assignments.tsv", 
                               sep='\t', fill=T, stringsAsFactors = F, header=T)


tax_paratha.V4V5 <- parathaData.V4V5 %>%
  dplyr::select(query.name, Kingdom, Phylum, Class, Order, Family, Genus, Species) 

parsed = seqinr::read.fasta(opts$fastaV4V5, #file('./input/SRR3225703.fasta'), 
                            as.string = TRUE,
                            forceDNAtolower = FALSE, whole.header = FALSE)

table = data.frame("ASV"=unlist(parsed), "query.name" = sapply(parsed, attr, 'name'), row.names=NULL)


tax_paratha.V4V5 <- merge(tax_paratha.V4V5, table)  %>%
  dplyr::group_by(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  dplyr::count(ASV, name="SRR3225703.V4V5") 

tax_paratha.V4V5 <- refine_species_names(tax_paratha.V4V5)

taxmatV4V5 <- tax_paratha.V4V5 %>% as.matrix
rownames(taxmatV4V5) <- paste0("x", tax_paratha.V4V5$ASV )
taxmatV4V5 <- taxmatV4V5[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]  
TAX_parathaa.V4V5 <- tax_table(taxmatV4V5)

otumatV4V5 <- as.matrix(tax_paratha.V4V5$SRR3225703.V4V5)
rownames(otumatV4V5) <- paste0("x", tax_paratha.V4V5$ASV)
colnames(otumatV4V5) <- "SRR3225703b.V4V5"
OTU_parathaa.V4V5 <- otu_table(otumatV4V5, taxa_are_rows = TRUE)

samp_parathaa.V4V5 <- data.frame(colnames(OTU_parathaa.V4V5), 
                                 rep("Parathaa", length(colnames(OTU_parathaa.V4V5))),
                                 rep("V4V5", length(colnames(OTU_parathaa.V4V5))))
colnames(samp_parathaa.V4V5) <- c("sampleID", "Taxonomy_type", "Region")
rownames(samp_parathaa.V4V5) <- samp_parathaa.V4V5$sampleID
str(samp_parathaa.V4V5)
SAMP_parathaa.V4V5 <- sample_data(samp_parathaa.V4V5)

ps1_parathaa.V4V5 <- phyloseq(OTU_parathaa.V4V5, TAX_parathaa.V4V5, SAMP_parathaa.V4V5)

## V1V2 -- Note that these are run on a "unique" seqs file, so the counts table needs to be imported separately

parathaData.V1V2 <- read.delim(
  opts$paraAssignV1V2, #"./output/20231205_MockV1V2/taxonomic_assignments.tsv", 
                               sep='\t', fill=T, stringsAsFactors = F, header=T)
v1v2.counts <- read.table(opts$V1V2Counts, #"input/SRR3225701.count_table", 
                          header = T)

tax_paratha.V1V2 <- parathaData.V1V2 %>%
  dplyr::select(query.name, Kingdom, Phylum, Class, Order, Family, Genus, Species) 



#Add in count table from dedup
tax_paratha.V1V2 <- left_join(tax_paratha.V1V2, v1v2.counts, by=c("query.name"="Representative_Sequence"))

parsed = seqinr::read.fasta(opts$fastaV1V2, as.string = TRUE,
                            forceDNAtolower = FALSE, whole.header = FALSE)

table = data.frame("ASV"=unlist(parsed), "query.name" = sapply(parsed, attr, 'name'), row.names=NULL)
#View(tax_paratha)

tax_paratha.V1V2 <- left_join(tax_paratha.V1V2, table) %>% 
  rename("SRR3225701.V1V2"=total)

#tax_paratha.V1V2 <- refine_species_names(tax_paratha.V1V2)
## "refine_species_names" takes too long on this, so going to manually change C. beijerinckii:

Cber <- tax_paratha.V1V2 %>% filter(Species == "Clostridium beijerinckii") %>% pull(query.name)
tax_paratha.V1V2[which(tax_paratha.V1V2$query.name %in% Cber), "Species"] <- "Clostridium beijerinckii/diolis"

taxmatV1V2 <- tax_paratha.V1V2 %>% as.matrix
rownames(taxmatV1V2) <- paste0("x", tax_paratha.V1V2$ASV )
taxmatV1V2 <- taxmatV1V2[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]  
TAX_parathaa.V1V2 <- tax_table(taxmatV1V2)

otumatV1V2 <- as.matrix(tax_paratha.V1V2$SRR3225701.V1V2)
rownames(otumatV1V2) <- paste0("x", tax_paratha.V1V2$ASV )
colnames(otumatV1V2) <- "SRR3225701b.V1V2"
OTU_parathaa.V1V2 <- otu_table(otumatV1V2, taxa_are_rows = TRUE)

samp_parathaa.V1V2 <- data.frame(colnames(OTU_parathaa.V1V2), 
                                 rep("Parathaa", length(colnames(OTU_parathaa.V1V2))),
                                 rep("V1V2", length(colnames(OTU_parathaa.V1V2))))
colnames(samp_parathaa.V1V2) <- c("sampleID", "Taxonomy_type", "Region")
rownames(samp_parathaa.V1V2) <- samp_parathaa.V1V2$sampleID
SAMP_parathaa.V1V2 <- sample_data(samp_parathaa.V1V2)

ps1_parathaa.V1V2 <- phyloseq(OTU_parathaa.V1V2, TAX_parathaa.V1V2, SAMP_parathaa.V1V2)


#this will be set to false... eventually this should be removed to clean up the code since we don't include this
# anaylsis in the final paper...
if(includeSpingo){
  ## Read in SPINGO results
  spingo.V4V5 <- read.delim("./output/20221101_MiSeqV4V5Mock/SRR3225703.summary.txt",
                            header = F, col.names = c("Species", "SRR3225703c"))
  spingo.V4V5$Species <- str_replace_all(spingo.V4V5$Species, "\\(", "")
  spingo.V4V5$SRR3225703c <- as.numeric(str_replace_all(spingo.V4V5$SRR3225703c, "\\)", ""))
  
  spingo.V4V5$Species <- str_replace_all(spingo.V4V5$Species, "_", " ")
  spingo.V4V5 <- refine_species_names(spingo.V4V5)
  
  taxmat <- as.matrix(spingo.V4V5$Species)
  rownames(taxmat) <- spingo.V4V5$Species 
  colnames(taxmat) <- "Species"
  TAX_spingo <- tax_table(taxmat)
  
  otumat <- as.matrix(spingo.V4V5$SRR3225703c)
  rownames(otumat) <- spingo.V4V5$Species
  colnames(otumat) <- "SRR3225703c"
  OTU_spingo <- otu_table(otumat, taxa_are_rows = TRUE)
  
  samp_spingo <- data.frame(colnames(OTU_spingo), rep("SPINGO", length(colnames(OTU_spingo))),
                            rep("V4V5", length(colnames(OTU_spingo))))
  colnames(samp_spingo) <- c("sampleID", "Taxonomy_type", "Region")
  rownames(samp_spingo) <- samp_spingo$sampleID
  str(samp_spingo)
  SAMP_spingo <- sample_data(samp_spingo)
  
  ps1_spingo.V4V5 <- phyloseq(OTU_spingo, TAX_spingo, SAMP_spingo)
  
  spingo.V1V2 <- read.delim("./output/20221219_MiSeqV1V2Mock/SRR3225701.summary.txt",
                            header = F, col.names = c("Species", "SRR3225701c"))
  spingo.V1V2$Species <- str_replace_all(spingo.V1V2$Species, "\\(", "")
  spingo.V1V2$SRR3225701c <- as.numeric(str_replace_all(spingo.V1V2$SRR3225701c, "\\)", ""))
  
  spingo.V1V2$Species <- str_replace_all(spingo.V1V2$Species, "_", " ")
  spingo.V1V2 <- refine_species_names(spingo.V1V2)
  
  taxmat <- as.matrix(spingo.V1V2$Species)
  rownames(taxmat) <- spingo.V1V2$Species 
  colnames(taxmat) <- "Species"
  TAX_spingo <- tax_table(taxmat)
  
  otumat <- as.matrix(spingo.V1V2$SRR3225701c)
  rownames(otumat) <- spingo.V1V2$Species
  colnames(otumat) <- "SRR3225701c"
  OTU_spingo <- otu_table(otumat, taxa_are_rows = TRUE)
  
  samp_spingo <- data.frame(colnames(OTU_spingo), rep("SPINGO", length(colnames(OTU_spingo))),
                            rep("V1V2", length(colnames(OTU_spingo))))
  colnames(samp_spingo) <- c("sampleID", "Taxonomy_type", "Region")
  rownames(samp_spingo) <- samp_spingo$sampleID
  str(samp_spingo)
  SAMP_spingo <- sample_data(samp_spingo)
  
  ps1_spingo.V1V2 <- phyloseq(OTU_spingo, TAX_spingo, SAMP_spingo)
  
  #load("output/mock_analysis/ps1_dada.V1V2.RData") # precalculated
  
  
  
  ps1_all<- merge_phyloseq(ps1_dada.V1V2, ps1_dada.V4V5, 
                           ps1_parathaa.V1V2, ps1_parathaa.V4V5, 
                           ps1_spingo.V1V2, ps1_spingo.V4V5)
}

if(includeSpingo==F){
  ps1_all<- merge_phyloseq(ps1_dada.V1V2, ps1_dada.V4V5, 
                           ps1_parathaa.V1V2, ps1_parathaa.V4V5)
}

ps1.com <- ps1_all
#saveRDS(pst.com, "ps1.com.RDS")

# We need to set Palette
taxic <- as.data.frame(ps1.com@tax_table) # this will help in setting large color options

#getPalette = colorRampPalette(brewer.pal(colourCount, ))  # change the palette as well as the number of colors will change according to palette.

taxic$OTU <- rownames(taxic) # Add the OTU ids from OTU table into the taxa table at the end.
colnames(taxic) # You can see that we now have extra taxonomy levels.

taxmat <- as.matrix(taxic) # convert it into a matrix.
new.tax <- tax_table(taxmat) # convert into phyloseq compatible file.
tax_table(ps1.com) <- new.tax # incroporate into phyloseq Object

# now edit the unclassified taxa
tax_table(ps1.com)[is.na(tax_table(ps1.com)[, "Phylum"]), "Phylum"] <- "p__"
tax_table(ps1.com)[is.na(tax_table(ps1.com)[, "Class"]), "Class"] <- "c__"
tax_table(ps1.com)[is.na(tax_table(ps1.com)[, "Order"]), "Order"] <- "o__"
tax_table(ps1.com)[is.na(tax_table(ps1.com)[, "Family"]), "Family"] <- "f__"
tax_table(ps1.com)[is.na(tax_table(ps1.com)[, "Genus"]), "Genus"] <- "Unknown"
tax_table(ps1.com)[which(tax_table(ps1.com)[, "Genus"]==""), "Genus"] <- "Unknown"
tax_table(ps1.com)[is.na(tax_table(ps1.com)[, "Species"]), "Species"] <- "Unknown"
tax_table(ps1.com)[which(tax_table(ps1.com)[, "Species"]=="UNCLASSIFIED"), "Species"] <- "Unknown"
tax_table(ps1.com)[which(tax_table(ps1.com)[, "Species"]=="AMBIGUOUS"), "Species"] <- "Unknown"


# it would be nice to have the Taxonomic names in italics.
# for that we set this

guide_italics <- guides(fill = guide_legend(label.theme = element_text(
  size = 15,
  face = "italic", colour = "Black", angle = 0
)))


## Now we need to plot at family level, we can do it as follows:

# first remove the phy_tree

ps1.com@phy_tree <- NULL

plotList <- list()

# Taxonomy plots
for(level in c("Genus","Species")){

  ps1.com.rel <- microbiome::transform(ps1.com, "compositional")
  ps1.com.rel.lev <- aggregate_rare(ps1.com.rel, level, detection = .0001/100, prevalence = 10/100)
  ps1.com.rel.lev.agg <- aggregate_top_taxa2(ps1.com.rel.lev, 14, level)

  if(level=="Genus")
    ps1.com.rel.lev.agg <- subset_samples(ps1.com.rel.lev.agg, Taxonomy_type!="SPINGO")
  
  #ps1.com.lev <- ps1.com

  taxa_names(ps1.com.rel.lev.agg) <- gsub("phBC6A52;", "phBC6A52;\n", taxa_names(ps1.com.rel.lev.agg))
  plot.composition.relAbun <- plot_composition(ps1.com.rel.lev.agg,
                                               #sample.sort="Taxonomy_type",
                                               ##sample.sort = "sampleID", 
                                               sample.sort = "Region",# for mock
                                               otu.sort = NULL,
                                               x.label = "sampleID") 
  #xlabs <- rep("", nrow(sample_data(ps1.com.rel.lev.agg)))
  #xlabs[ceiling(nrow(sample_data(ps1.com.rel.lev.agg))/4)] <- "DADA2"
  #xlabs[ceiling(nrow(sample_data(ps1.com.rel.lev.agg))*3/4)+1] <- "parathaa"
  if(level=="Genus" | includeSpingo==F)
    xlabs <- c("DADA2 V1V2", "Parathaa V1V2", "DADA2 V4V5", "Parathaa V4V5")
  if(level=="Species" & includeSpingo==T)
    xlabs <- c("DADA2 V1V2", "Parathaa V1V2", "SPINGO V1V2", "DADA2 V4V5", "Parathaa V4V5", "SPINGO V4V5") # for mock
  plot.composition.relAbun <- plot.composition.relAbun + theme(legend.position = "bottom") 
  plot.composition.relAbun$data$Tax <- relevel(plot.composition.relAbun$data$Tax, "Other")
  plot.composition.relAbun$data$Tax <- relevel(plot.composition.relAbun$data$Tax, "Unknown")
  cols <- c("white", "grey", stepped2(n=length(levels(plot.composition.relAbun$data$Tax))-2))
  plot.composition.relAbun <- plot.composition.relAbun +  theme_bw()  +  scale_fill_manual(level, values = cols)#scale_fill_brewer(palette="Paired")  
  #plot.composition.relAbun <- plot.composition.relAbun + theme(axis.text.x = element_text(angle = 90)) 
  plot.composition.relAbun <- plot.composition.relAbun + ylab("Relative Abundance") +
    guides(fill= 
             guide_legend(
               ncol=1,
               label.theme = element_text(face="italic", size=14),
               title.theme = element_text(size=18),
               )
           )
  plot.composition.relAbun <- plot.composition.relAbun +   theme(axis.text.x = element_text(size=12, angle=45, hjust=1)) + scale_x_discrete(labels=xlabs) 

  if(level!="Species")
    ggsave(filename=paste0(opts$o, "/Figures/", level, "mock_taxonomy.png"), plot.composition.relAbun, width=6, height=6)
  if(level=="Species")
    ggsave(filename=paste0(opts$o, "/Figures/", level, "mock_taxonomy.png"), plot.composition.relAbun, width=10, height=6)
  
  plotList[[level]] <- plot.composition.relAbun
}

## Panel Plot ##
p1 <- ggarrange(plotList[["Genus"]], plotList[["Species"]] , ncol=2, labels = c("B", "C"), align="h", widths=c(1,1.5)) 
ggsave(filename = paste0(opts$o, "/Figures/", "mock_genus_species_panel.png"), plot= p1, height = 6, width=14, dpi=600)
## Panel Plot ##
p1 <- ggarrange(plotList[["Species"]], plotList[["Genus"]], ncol=2, labels = c("B", "C"), align="h", widths=c(1.5,1)) 
ggsave(filename = paste0(opts$o, "/Figures/", "mock_species_genus_panel.png"), plot= p1, height = 6, width=14, dpi=600)


## Heatmap of mock community
level <- "Species"
ps1.com.rel <- microbiome::transform(ps1.com, "compositional")
ps1.com.rel.lev <- aggregate_rare(ps1.com.rel, level, detection = .1/100, prevalence = 10/100)
ps1.com.rel.lev.agg <- aggregate_taxa(ps1.com.rel.lev, level)
forHeat <- t(as.matrix(otu_table(ps1.com.rel.lev.agg)))

if(includeSpingo){
  rownames(forHeat) <- c("V1V2 DADA2", "V4V5 DADA2", 
                         "V1V2 Parathaa", "V4V5 Parathaa", 
                         "V1V2 SPINGO", "V4V5 SPINGO")
  forHeat <- forHeat[c(1,3,5,2,4,6),]
}

if(includeSpingo==F){
  rownames(forHeat) <- c("V1V2 DADA2", "V4V5 DADA2", 
                         "V1V2 Parathaa", "V4V5 Parathaa")
  forHeat <- forHeat[c(1,3,2,4),]
}

included <- c("Actinomyces odontolyticus",
              "Bacillus cereus",
              "Bacillus cereus;Bacillus thuringiensis",
              "Bacillus anthracis;Bacillus cereus;Bacillus thuringiensis",
              "Bacteroides vulgatus"                                     ,
              "Clostridium diolis"                                        ,
              "Clostridium sensu stricto 1 diolis"                        ,
              "Clostridium beijerinckii",
              "Clostridium beijerinckii/diolis",
              "Clostridium beijerinckii;Clostridium diolis",
              "Cutibacterium acnes",
              "Deinococcus radiodurans"                                  ,
              "Enterococcus canis;Enterococcus faecalis"         ,
              "Enterococcus canis;Enterococcus casseliflavus;Enterococcus faecalis;Enterococcus gallinarum",
              "Enterococcus faecalis"                                    ,
              "Lactobacillus gasseri;Lactobacillus johnsonii"            ,
              "Lactobacillus gasseri",
              "Listeria innocua;Listeria monocytogenes"                   ,
              "Listeria monocytogenes"                                   ,
              "Neisseria meningitidis"                                   ,
              "Propionibacterium acnes", 
              "Pseudomonas aeruginosa"                                  ,
              "Staphylococcus aureus"                                    ,
              "Staphylococcus epidermidis"  ,
              "Staphylococcus aureus;Staphylococcus epidermidis"         ,  
              "Streptococcus agalactiae"                                 ,
              "Streptococcus agalactiae;Streptococcus suis",  
              "Streptococcus mutans"                                     ,
              "Streptococcus pneumoniae",
              "Enterococcus canis;Enterococcus casseliflavus;Enterococcus faecalis;Enterococcus gallinarum;Enterococcus Unclassified;Melissococcus plutonius",
              "Bacillus anthracis;Bacillus cereus;Bacillus phage phBC6A52;Bacillus thuringiensis;Bacillus Unclassified",
              "Bacillus cereus;Bacillus thuringiensis;Bacillus Unclassified;Bacillus wiedmannii",
              "Streptococcus agalactiae;Streptococcus phage 10750.4")
others <- c("Ambiguous", "Unknown", "Other")
included.df <- data.frame(rep("Not included", length(colnames(forHeat))))
colnames(included.df) <- "Mock Community"
included.df[which(colnames(forHeat) %in% included), "Mock Community"] <- "Included"
included.df[which(colnames(forHeat) %in% others), "Mock Community"] <- "Other"
rownames(included.df) <- colnames(forHeat)

type <- gsub("s\\d+_", "", colnames(forHeat))
col_fun = viridis(3)
ha = columnAnnotation(
  df = included.df[,1], show_annotation_name=FALSE, 
  annotation_legend_param = list(title = "Mock Community"),
  col=list(df = c("Included" =  viridis(4)[2],"Not included" = viridis(4)[3], "Other" = viridis(4)[1]))
)
region.df <- data.frame(c(rep("V1V2", length(rownames(forHeat))/2), rep("V4V5", length(rownames(forHeat))/2)))
colnames(region.df) <- "16S Region"
rownames(region.df) <- rownames(forHeat)
  
ra = rowAnnotation(
  df = region.df[,1], show_annotation_name=FALSE, 
  annotation_legend_param = list(title = "16S Region"),
  col=list(df = c("V1V2" =  turbo(7)[2],"V4V5" = turbo(7)[4]))
)
col_fun <- c("grey", "grey", "grey",   rev(magma(98)))

### Anything < 0.1% set to 0
forHeat <- forHeat*(forHeat>0.001)
###manually fix heatmap column names...

colnames(forHeat)[1] <- "Bacillus anthracis;B. cereus;B. phage phBC6A52;B. thuringiensis;B. Unclassified"
colnames(forHeat)[2] <- "Bacillus anthracis;B. cereus;B. thuringiensis" 
colnames(forHeat)[4] <- "Bacillus cereus;B. thuringiensis;B. Unclassified;B. wiedmannii"
colnames(forHeat)[7] <- "Enterococcus canis;E. casseliflavus;E. faecalis;E. gallinarum;E. Unclassified;Melissococcus plutonius"
colnames(forHeat)[9] <- "Lactobacillus gasseri;L. johnsonii"
colnames(forHeat)[10] <- "Listeria innocua;L. monocytogenes"
colnames(forHeat)[15] <- "Staphylococcus aureus;S. epidermidis"
colnames(forHeat)[18] <- "Streptococcus agalactiae;S. phage 10750.4"

png(paste0(opts$o, "/Fig_2A_MockHeatmap.png"),width=9.5,height=4.5,units="in", res=600)
ht1 <- Heatmap(sqrt(forHeat), name = "sqrt(Rel.\nabundance)", col=col_fun,
               column_names_max_height = unit(15, "cm"),
               cluster_columns = FALSE, cluster_rows=FALSE, column_names_rot = 45, 
               column_names_side = "top", column_names_gp = grid::gpar(fontsize = 7), 
               column_split = included.df[,"Mock Community"], column_title=NULL, top_annotation = ha,
               row_split = region.df[,"16S Region"], row_title=NULL, right_annotation = ra,
               row_names_gp = grid::gpar(fontsize = 10))  %v% NULL

draw(ht1)

dev.off()

