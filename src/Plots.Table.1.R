## Compare taxonomies from DADA2 vs PARATHAA


require(docopt)

'Usage:
  Plots.Table.1.R [-p <parathaa_PATH> --dada_db <dada2_db> --dada_db_sp <dada2_db_sp> -t <input_taxonomy> -o <output> --paraAssignV4V5 <V4V5_taxonomy_file> --paraAssignV1V2 <V1V2_taxonomy_file> --queryV4V5 <V4V5 seqs> --queryV1V2 <V1V2 seqs>]
  
  
Options:
  -p directory where parathaa github repo is cloned
  --dada_db location of DADA2 db
  --dada_db_sp location of DADA2 species db
  --paraAssignV4V5 parathaa assignments for V4V5
  --paraAssignV1V2 parathaa assignments for V1V2
  --queryV4V5 query sequence file for V4V5
  --queryV1V2 query sequence file for V1V2
  -t input taxonomy
  -o output
  ]' -> doc


opts <- docopt(doc)

parathaaDir <- (opts$p)

library(phyloseq) 
library(dplyr) 
library(stringr)
library(ggtree)
library(treeio)
library(dada2)
library(ggplot2)
library(castor)
library(phytools)
#library(flextable)
suppressPackageStartupMessages(library(seqinr))
source(file.path(parathaaDir, "parathaa/utility/SILVA.species.editor.dev.R"))

##this should be fine I think because we won't be installing the manuscript tools.
source("src/performance.table.R")



##################################
## Synthetic community analysis ##
##################################

## Specify dada2 dbs

#DADAdb <- "input/20230111.silva.seed_v138_1.ng.dada.fasta"
#DADAdb.sp <- "input/20231130_silva.seed_v138_1.ng.dada.sp.fasta"
DADAdb <- opts$dada_db
DADAdb.sp <- opts$dada_db_sp



## Define Function ##
run.synthetic.data <- function(parathaaFile, sequenceFile, regionName, outputDir, dadaAllowMult = FALSE){
  
  dir.create(outputDir, recursive = T, showWarnings = F)
  
  # Read in parathaa data
  parathaaData <- read.delim(
    parathaaFile,
    sep='\t', fill=T, stringsAsFactors = F, header=T)
  tax_parathaa <- parathaaData %>%
    dplyr::select(query.name, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    group_by(query.name) 
  
  taxmat <- tax_parathaa %>% as.matrix
  rownames(taxmat) <- tax_parathaa$query.name
  taxmat <- taxmat[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]  
  TAX_parathaa <- tax_table(taxmat)
  
  otutab <- tax_parathaa %>% 
    dplyr::group_by(query.name, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
    dplyr::count(query.name, name="Parathaa") 
  otumat <- as.matrix(otutab$Parathaa)
  rownames(otumat) <- otutab$query.name
  colnames(otumat) <- paste0("Parathaa.", regionName)
  OTU_parathaa <- otu_table(otumat, taxa_are_rows = TRUE)
  
  samp_parathaa <- data.frame(colnames(OTU_parathaa), rep("parathaa", length(colnames(OTU_parathaa))),
                              rep(regionName, length(colnames(OTU_parathaa))))
  colnames(samp_parathaa) <- c("sampleID", "Taxonomy_type", "Region")
  rownames(samp_parathaa) <- samp_parathaa$sampleID
  str(samp_parathaa)
  SAMP_parathaa <- sample_data(samp_parathaa)
  
  ps1_parathaa <- phyloseq(OTU_parathaa, TAX_parathaa, SAMP_parathaa)
  print(ps1_parathaa)
  
  ## Assign taxonomy with DADA2
  ## First, get names and sequences from fasta file
  getNames <- read.fasta(file(sequenceFile), as.string = TRUE,
                         forceDNAtolower = FALSE, whole.header = FALSE)
  names1 <- str_split(getName(getNames), "\t", simplify=TRUE)
  names1 <- names1[,1] %>%
    str_remove(">")
  name.df <- data.frame("sequence" = unlist(getSequence(getNames, as.string=T)), taxaIDs = names1)
  
  nChars2 <- name.df %>% filter(str_detect(sequence, "N")) %>% pull(taxaIDs)
  
  ## Next, assign taxonomy to genus level with DADA2 (takes a few minutes)
  set.seed(3874)
  taxa <- assignTaxonomy(sequenceFile, 
                         DADAdb,
                         multithread=TRUE)
  
  
  ## Remove sequences with undefined ("N") bases, store until after species assignment
  taxa.test <- as.data.frame(taxa)
  taxa.test$taxaIDs <- names1
  nChars <- grep("N", rownames(taxa.test))
  print(paste("Removing", length(nChars), "sequences with N bases"))
  withNbases <- taxa.test[nChars,]
  taxa <- taxa[-nChars,]
  
  ## Perform species assignment with DADA2 (takes a few minutes)
  taxa.sp <- addSpecies(taxa, DADAdb.sp, allowMultiple = dadaAllowMult)
  
  #taxa.sp.mult <- addSpecies(taxa, DADAdb.sp, allowMultiple = T)
  
  ## Add in reference IDs and taxonomy from sequences with "N" bases
  tax_dada <- as.data.frame(taxa.sp)
  tax_dada$sequence <- str_split(rownames(tax_dada), "\\.", simplify=TRUE)[,1]
  getnamSubset <-name.df %>% filter(sequence %in% tax_dada$sequence)
  tax_dada2 <- cbind(tax_dada, "taxaIDs" =getnamSubset$taxaIDs)
  tax_dada2 <- full_join(tax_dada2, withNbases)
  rownames(tax_dada2) <- tax_dada2[,"taxaIDs"]
  
  tax_dada3 <- tax_dada2 %>%
    select(Kingdom, Phylum, Class, Order, Family, Genus, Species) 
  
  #New Species creation:
  tax_dada3 <- tax_dada3 %>% 
    rowwise() %>% 
    mutate(Species = if_else(!is.na(Species), 
                             paste( paste(Genus), str_split(Species, "/",simplify = T), collapse =";"), 
                             NA)
    ) %>% 
    as.matrix()
  rownames(tax_dada3) <- tax_dada2[,"taxaIDs"]
  
  
  ## Place DADA2 taxonomies into phyloseq object
  TAX_dada <- tax_table(tax_dada3)
  
  otutab <- as.data.frame(tax_dada2) %>% 
    select(taxaIDs, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    dplyr::group_by(taxaIDs, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
    dplyr::count(taxaIDs, name="DADA2") 
  otumat <- as.matrix(otutab$DADA2)
  rownames(otumat) <- otutab$taxaIDs
  colnames(otumat) <- paste0("DADA2.", regionName)
  OTU_dada <- otu_table(otumat, taxa_are_rows = TRUE)
  
  samp_dada <- data.frame(colnames(OTU_dada), rep("DADA2", length(colnames(OTU_dada))),
                          rep(regionName, length(colnames(OTU_dada))))
  colnames(samp_dada) <- c("sampleID", "Taxonomy_type", "Region")
  rownames(samp_dada) <- samp_dada$sampleID
  str(samp_dada)
  SAMP_dada <- sample_data(samp_dada)
  
  ps1_dada <- phyloseq(OTU_dada, TAX_dada, SAMP_dada)
  
  print(ps1_dada)
  
  
  ############################################
  ### Assess performance on Synthetic Data ###
  ############################################
  
  #Get taxdata from SILVA:

  inFileTaxdata <- opts$t
  
  taxdata <- read.table(inFileTaxdata , header=T, fill=TRUE,sep='\t', quote="")
  taxdata <- taxdata %>%
    unite("AccID", c("primaryAccession", "start", "stop"), sep=".", remove=F)
  taxdata <- taxdata %>%
    mutate(taxonomy=paste0(path, organism_name))
  taxdata <- SILVA.species.editor.DADA(taxdata, "taxonomy")
  
  taxdata <- taxdata %>%
    select(AccID, primaryAccession, start, stop, taxonomy) %>%
    separate(col=taxonomy, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
  
  ## Additional changes (getting rid of subspecies)
  taxdata <- SILVA.species.editor(taxdata)
  taxdata <- taxdata %>% filter(!is.na(Species))
  
  
  ## Read in seed db for exclusion
  inFileSeedDB <- "input/silva.seed_v138_1.tax"
  SeedTax <- read.table(inFileSeedDB , header=F, fill=TRUE,sep='\t')
  
  SeedTax <- SeedTax %>%
    separate(col=V1, into=c("primaryAccession", "ArbID"), sep="\\.") %>%
    separate(col=V2, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep=";") %>%
    filter(Kingdom=="Bacteria" & !is.na(Genus) & Genus!="") 
  
  taxdata_seedless <- taxdata %>% 
    filter(!primaryAccession %in% SeedTax$primaryAccession) %>%
    ## Subset to Genera in Seed for now
    filter(Genus %in% unique(SeedTax$Genus))
  
  ## Make synthetic comparison dataset
  synth.parathaa<- as.data.frame(tax_table(ps1_parathaa))
  synth.parathaa$AccID <- rownames(synth.parathaa)
  
  synth.parathaa2 <- left_join(synth.parathaa, taxdata, by="AccID")
  synth.parathaa2 <- synth.parathaa2 %>% 
    mutate(Species.x = unlist(lapply(str_split(Species.x, ";"), FUN=function(x) paste0(word(x,1,2), collapse = ";" ))))
  synth.parathaa2 <- synth.parathaa2 %>% 
    mutate(Species.x = ifelse(Species.x=="NA", NA, Species.x),
           Genus.x = ifelse(Genus.x=="", NA, Genus.x))
  synth.parathaa3 <- synth.parathaa2 %>% 
    dplyr::rowwise() %>%
    mutate(Flag = ifelse(is.na(Species.x), NA, word(Species.y, 1, 2) %in% str_split(Species.x, ";", simplify = T)),
           Flag.genus = ifelse(is.na(Genus.x), NA, Genus.y %in% str_split(Genus.x, ";", simplify = T))
    )
  
  synth.dada <- as.data.frame(tax_table(ps1_dada))
  synth.dada$AccID <- rownames(synth.dada)
  
  synth.dada2 <- left_join(synth.dada, taxdata, by="AccID")
  synth.dada2 <- synth.dada2 %>% 
    mutate(Species.x = unlist(lapply(str_split(Species.x, ";"), FUN=function(x) paste0(word(x,1,2), collapse = ";" ))))
  synth.dada2 <- synth.dada2 %>% 
    mutate(Species.x = ifelse(Species.x=="NA", NA, Species.x),
           Genus.x = ifelse(Genus.x=="", NA, Genus.x))
  synth.dada3 <- synth.dada2 %>% 
    dplyr::rowwise() %>%
    mutate(Flag = ifelse(is.na(Species.x), NA, word(Species.y, 1, 2) %in% str_split(Species.x, ";", simplify = T)),
           Flag.genus = ifelse(is.na(Genus.x), NA, Genus.y %in% str_split(Genus.x, ";", simplify = T))
    )
  
  
  compare.synth <- dplyr::full_join(synth.dada3, synth.parathaa3, by="AccID")
  compare.synth <- compare.synth %>% 
    mutate(Species.silva = word(Species.y.y, 1, 2)) %>%
    rename(Species.parathaa = Species.x.y,
           Species.dada = Species.x.x,
           Genus.dada = Genus.x.x,
           Genus.parathaa = Genus.x.y,
           Genus.silva = Genus.y.y)
  
 ### Remove seqs with N characters: 
  compare.synth <- compare.synth %>% filter(!AccID %in% nChars2)
  save(compare.synth, file= file.path(outputDir, paste0(regionName, "_full_comparisons.RData")))
  
  t1 <- performance.table(compare.synth, "Species")
  write.table(t1, file= file.path(outputDir, paste0(regionName, "_Species_performance.tsv")), sep="\t", col.names = NA)
  t2 <- performance.table(compare.synth, "Genus")
  write.table(t2, file= file.path(outputDir, paste0(regionName, "_Genus_performance.tsv")), sep="\t", col.names = NA)
  
  
  compare.genus <- data.frame()
  ## Both uniquely correct:
  compare.genus["Dada2_correct", "Parathaa_correct"] <- nrow(compare.synth %>% filter(Genus.dada==Genus.silva & Genus.parathaa==Genus.silva) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
  ## Paratha uniquely correct, dada2 incorrect:
  compare.genus["Dada2_incorrect", "Parathaa_correct"] <- nrow(compare.synth %>% filter(!Flag.genus.x & !is.na(Genus.dada)  & Genus.parathaa==Genus.silva) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
  ## Paratha uniquely correct, dada2 unassigned:
  compare.genus["Dada2_unassigned", "Parathaa_correct"] <- nrow(compare.synth %>% filter(is.na(Flag.genus.x)  & Genus.parathaa==Genus.silva) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
  
  ##Paratha partly correct, dada2 uniquely correct:
  compare.genus["Dada2_correct", "Parathaa_1toMany"] <- nrow(compare.synth %>% filter(Genus.dada==Genus.silva & Flag.genus.y & Genus.parathaa!=Genus.silva) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
  ##Paratha partly correct, dada2 incorrect:
  compare.genus["Dada2_incorrect", "Parathaa_1toMany"] <- nrow(compare.synth %>% filter(!Flag.genus.x & Flag.genus.y & Genus.parathaa!=Genus.silva) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
  ##Paratha partly correct, dada2 unassigned:
  compare.genus["Dada2_unassigned", "Parathaa_1toMany"] <- nrow(compare.synth %>% filter(is.na(Flag.genus.x) & Flag.genus.y & Genus.parathaa!=Genus.silva) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
  
  ##Paratha incorrect, dada2 uniquely correct:
  compare.genus["Dada2_correct", "Parathaa_incorrect"] <- nrow(compare.synth %>% filter(Genus.dada==Genus.silva & !Flag.genus.y ) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
  ##Both incorrect:
  compare.genus["Dada2_incorrect", "Parathaa_incorrect"] <- nrow(compare.synth %>% filter(!Flag.genus.x & !Flag.genus.y ) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
  ##Paratha incorrect, dada2 unassigned:
  compare.genus["Dada2_unassigned", "Parathaa_incorrect"] <- nrow(compare.synth %>% filter(is.na(Flag.genus.x) & !Flag.genus.y ) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
  
  ##Paratha unassigned, dada2 uniquely correct:
  compare.genus["Dada2_correct", "Parathaa_unassigned"] <- nrow(compare.synth %>% filter(Genus.dada==Genus.silva & is.na(Flag.genus.y) ) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
  ##Paratha unassigned, dada2 incorrect:
  compare.genus["Dada2_incorrect", "Parathaa_unassigned"] <- nrow(compare.synth %>% filter(!Flag.genus.x & is.na(Flag.genus.y)) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
  ##Both unassigned:
  compare.genus["Dada2_unassigned", "Parathaa_unassigned"] <- nrow(compare.synth %>% filter(is.na(Flag.genus.x) & is.na(Flag.genus.y) ) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
  
  sum(compare.genus)
  compare.genus$Sum <- apply(compare.genus, 1, sum)
  compare.genus["Sum",] <- apply(compare.genus, 2, sum)
  
  compare.genus[,"Pct"] <- round(compare.genus[,"Sum"] / compare.genus["Sum", "Sum"], 3)
  compare.genus["Pct",] <- round(compare.genus["Sum",] / compare.genus["Sum", "Sum"], 3)
  
  cat("Genus comparison:", sep = "\n", file = file.path(outputDir, paste0(regionName, "synthetic_output.txt")), append = F)
  capture.output(t(compare.genus), file = file.path(outputDir, paste0(regionName, "synthetic_output.txt")), append = T)
  
  ###################################
  #### Species-level Summary ########
  ###################################

  
  ## Distinguish between query species "included" in the database and those "novel" to the database
  taxdata_seed <- taxdata %>% 
    filter(primaryAccession %in% SeedTax$primaryAccession)
  taxdata_SP <- word(taxdata_seed$Species, 1, 2)
  
  compare.synth.full <- compare.synth
  compare.synth.novel <- compare.synth.full %>% filter(!Species.silva %in% taxdata_SP)
  compare.synth.included <- compare.synth.full %>% filter(Species.silva %in% taxdata_SP)
  
  sp.list <- list(compare.synth.full, compare.synth.novel, compare.synth.included)
  names(sp.list) <- c("All queries", "Novel Queries Only", "Included Queries Only")
  
  for(ind in 1:length(sp.list)){
    comp1 <- sp.list[[ind]]
    compare.species <- data.frame()
    
    if(dadaAllowMult==F){
      ## Both uniquely correct:
      compare.species["Dada2_correct", "Parathaa_correct"] <- nrow(comp1 %>% filter(Species.dada==Species.silva & Species.parathaa==Species.silva) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
      ## Paratha uniquely correct, dada2 incorrect:
      compare.species["Dada2_incorrect", "Parathaa_correct"] <- nrow(comp1 %>% filter(!Flag.x & !is.na(Species.dada)  & Species.parathaa==Species.silva) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
      ## Paratha uniquely correct, dada2 unassigned:
      compare.species["Dada2_unassigned", "Parathaa_correct"] <- nrow(comp1 %>% filter(is.na(Flag.x)  & Species.parathaa==Species.silva) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
      
      ##Paratha partly correct, dada2 uniquely correct:
      compare.species["Dada2_correct", "Parathaa_1toMany"] <- nrow(comp1 %>% filter(Species.dada==Species.silva & Flag.y & Species.parathaa!=Species.silva) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
      ##Paratha partly correct, dada2 incorrect:
      compare.species["Dada2_incorrect", "Parathaa_1toMany"] <- nrow(comp1 %>% filter(!Flag.x & Flag.y & Species.parathaa!=Species.silva) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
      ##Paratha partly correct, dada2 unassigned:
      compare.species["Dada2_unassigned", "Parathaa_1toMany"] <- nrow(comp1 %>% filter(is.na(Flag.x) & Flag.y & Species.parathaa!=Species.silva) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
      
      ##Paratha incorrect, dada2 uniquely correct:
      compare.species["Dada2_correct", "Parathaa_incorrect"] <- nrow(comp1 %>% filter(Species.dada==Species.silva & !Flag.y ) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
      ##Both incorrect:
      compare.species["Dada2_incorrect", "Parathaa_incorrect"] <- nrow(comp1 %>% filter(!Flag.x & !Flag.y ) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
      ##Paratha incorrect, dada2 unassigned:
      compare.species["Dada2_unassigned", "Parathaa_incorrect"] <- nrow(comp1 %>% filter(is.na(Flag.x) & !Flag.y ) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
      
      ##Paratha unassigned, dada2 uniquely correct:
      compare.species["Dada2_correct", "Parathaa_unassigned"] <- nrow(comp1 %>% filter(Species.dada==Species.silva & is.na(Flag.y) ) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
      ##Paratha unassigned, dada2 incorrect:
      compare.species["Dada2_incorrect", "Parathaa_unassigned"] <- nrow(comp1 %>% filter(!Flag.x & is.na(Flag.y)) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
      ##Both unassigned:
      compare.species["Dada2_unassigned", "Parathaa_unassigned"] <- nrow(comp1 %>% filter(is.na(Flag.x) & is.na(Flag.y) ) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
    }
    
    if(dadaAllowMult==T){
      ## Both uniquely correct:
      compare.species["Dada2_correct", "Parathaa_correct"] <- nrow(comp1 %>% filter(Species.dada==Species.silva & Species.parathaa==Species.silva) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
      ## Paratha uniquely correct, dada2 partly correct:
      compare.species["Dada2_1toMany", "Parathaa_correct"] <- nrow(comp1 %>% filter(Flag.x & Species.dada!=Species.silva  & Species.parathaa==Species.silva) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
      ## Paratha uniquely correct, dada2 incorrect:
      compare.species["Dada2_incorrect", "Parathaa_correct"] <- nrow(comp1 %>% filter(!Flag.x & !is.na(Species.dada)  & Species.parathaa==Species.silva) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
      ## Paratha uniquely correct, dada2 unassigned:
      compare.species["Dada2_unassigned", "Parathaa_correct"] <- nrow(comp1 %>% filter(is.na(Flag.x)  & Species.parathaa==Species.silva) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
      
      ##Paratha partly correct, dada2 uniquely correct:
      compare.species["Dada2_correct", "Parathaa_1toMany"] <- nrow(comp1 %>% filter(Species.dada==Species.silva & Flag.y & Species.parathaa!=Species.silva) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
      ## Paratha partly correct, dada2 partly correct:
      compare.species["Dada2_1toMany", "Parathaa_1toMany"] <- nrow(comp1 %>% filter(Flag.x & Species.dada!=Species.silva  & Flag.y & Species.parathaa!=Species.silva) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
      ##Paratha partly correct, dada2 incorrect:
      compare.species["Dada2_incorrect", "Parathaa_1toMany"] <- nrow(comp1 %>% filter(!Flag.x & Flag.y & Species.parathaa!=Species.silva) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
      ##Paratha partly correct, dada2 unassigned:
      compare.species["Dada2_unassigned", "Parathaa_1toMany"] <- nrow(comp1 %>% filter(is.na(Flag.x) & Flag.y & Species.parathaa!=Species.silva) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
      
      ##Paratha incorrect, dada2 uniquely correct:
      compare.species["Dada2_correct", "Parathaa_incorrect"] <- nrow(comp1 %>% filter(Species.dada==Species.silva & !Flag.y ) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
      ## Paratha incorrect, dada2 partly correct:
      compare.species["Dada2_1toMany", "Parathaa_incorrect"] <- nrow(comp1 %>% filter(Flag.x & Species.dada!=Species.silva  & !Flag.y ) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
      ##Both incorrect:
      compare.species["Dada2_incorrect", "Parathaa_incorrect"] <- nrow(comp1 %>% filter(!Flag.x & !Flag.y ) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
      ##Paratha incorrect, dada2 unassigned:
      compare.species["Dada2_unassigned", "Parathaa_incorrect"] <- nrow(comp1 %>% filter(is.na(Flag.x) & !Flag.y ) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
      
      ##Paratha unassigned, dada2 uniquely correct:
      compare.species["Dada2_correct", "Parathaa_unassigned"] <- nrow(comp1 %>% filter(Species.dada==Species.silva & is.na(Flag.y) ) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
      ## Paratha unassigned, dada2 partly correct:
      compare.species["Dada2_1toMany", "Parathaa_unassigned"] <- nrow(comp1 %>% filter(Flag.x & Species.dada!=Species.silva  & is.na(Flag.y) ) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
      ##Paratha unassigned, dada2 incorrect:
      compare.species["Dada2_incorrect", "Parathaa_unassigned"] <- nrow(comp1 %>% filter(!Flag.x & is.na(Flag.y)) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
      ##Both unassigned:
      compare.species["Dada2_unassigned", "Parathaa_unassigned"] <- nrow(comp1 %>% filter(is.na(Flag.x) & is.na(Flag.y) ) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
    }
    
    
    sum(compare.species)
    compare.species$Sum <- apply(compare.species, 1, sum)
    compare.species["Sum",] <- apply(compare.species, 2, sum)
    
    compare.species[,"Pct"] <- round(compare.species[,"Sum"] / compare.species["Sum", "Sum"], 3)
    compare.species["Pct",] <- round(compare.species["Sum",] / compare.species["Sum", "Sum"], 3)
    
    
    
    
    cat("Species Output:", sep = "\n", file = file.path(outputDir, paste0(regionName, "synthetic_output.txt")), append = T)
    cat(names(sp.list)[ind] , sep = "\n", file = file.path(outputDir, paste0(regionName, "synthetic_output.txt")), append = T)
    capture.output(t(compare.species),  file = file.path(outputDir, paste0(regionName, "synthetic_output.txt")), append = T)

    
  }
  
}

##### CALL FUNCTION ######
run.synthetic.data(parathaaFile = opts$paraAssignV4V5, 
                   sequenceFile = opts$queryV4V5,
                   regionName = "V4V5", 
                   outputDir=paste(opts$o, "/Figures/synth_mult_arc/", sep=""),
                   dadaAllowMult = T)

run.synthetic.data(parathaaFile = opts$paraAssignV1V2, 
                   sequenceFile = "input/SILVAsubsample_SeedGenera_V1V2.pcr.fasta", 
                   regionName = "V1V2", 
                   outputDir=paste(opts$o, "/Figures/synth_mult_arc", sep=""),
                   dadaAllowMult = T)

run.synthetic.data(parathaaFile = opts$paraAssignV4V5, 
                   sequenceFile = opts$queryV4V5,
                   regionName = "V4V5", 
                   outputDir=paste(opts$o, "/Figures/synth_nomult_arc/", sep=""),
                   dadaAllowMult = F)

run.synthetic.data(parathaaFile = opts$paraAssignV1V2, 
                   sequenceFile = "input/SILVAsubsample_SeedGenera_V1V2.pcr.fasta", 
                   regionName = "V1V2", 
                   outputDir=paste(opts$o, "/Figures/synth_nomult_arc", sep=""),
                   dadaAllowMult = F)

##### MAKE NICE-LOOKING OUTPUT TABLE #####
# Below code will make microsoft word tables. Currently not run in the manuscript workflow. 
# ### Dada2 without AllowMultiple ###
# V4V5.tab <- read.table("output/Figures/synth_nomult_arc/V4V5_Species_performance.tsv", row.names=1, header = T, sep='\t')
# V1V2.tab <- read.table("output/Figures/synth_nomult_arc/V1V2_Species_performance.tsv", row.names=1, header=T, sep='\t')
# 
# V4V5.genus.tab <- read.table("output/Figures/synth_nomult_arc/V4V5_Genus_performance.tsv", row.names=1, header = T, sep='\t')
# V1V2.genus.tab <- read.table("output/Figures/synth_nomult_arc/V1V2_Genus_performance.tsv", row.names=1, header=T, sep='\t')
# 
# 
# summary.tab1 <- cbind(rownames(V1V2.tab), V1V2.tab, V4V5.tab, V1V2.genus.tab, V4V5.genus.tab)
# colnames(summary.tab1) <- c("Metric", "Species.V1V2.Parathaa", "Species.V1V2.DADA2", 
#                             "Species.V4V5.Parathaa", "Species.V4V5.DADA2", 
#                             "Genus.V1V2.Parathaa", "Genus.V1V2.DADA2", 
#                             "Genus.V4V5.Parathaa", "Genus.V4V5.DADA2")
# 
# ft1 <- flextable(summary.tab1) |>
#   separate_header() 
# border <- fp_border_default()
# ft1 <- vline(ft1, j = c('Metric','Species.V4V5.DADA2'), border = border, part = "all")
# ft1<- hline(ft1, i=4)
# ft1
# save_as_docx(ft1,  path = "output/Figures/synth_nomult_arc/Synth_performance.docx")
# 
# 
# ### Dada2 with AllowMultiplt ###
# V4V5.multdada <- read.table("output/Figures/synth_mult_arc/V4V5_Species_performance.tsv", row.names=1, header=T, sep='\t')
# V1V2.multdada <- read.table("output/Figures/synth_mult_arc/V1V2_Species_performance.tsv", row.names=1, header=T, sep='\t')
# 
# summary.tab2 <- cbind(rownames(V1V2.multdada), V1V2.multdada, V4V5.multdada)
# colnames(summary.tab2) <- c("Metric", "V1V2.Parathaa", "V1V2.DADA2", "V4V5.Parathaa", "V4V5.DADA2")
# 
# ft2 <- flextable(summary.tab2) |>
#   separate_header() 
# border <- fp_border_default()
# ft2 <- vline(ft2, j = c('Metric','V1V2.DADA2'), border = border, part = "all")
# ft2<- hline(ft2, i=4)
# ft2
# 
# save_as_docx(ft2,  path = "output/Figures/synth_mult_arc/Synth_performance_mult.docx")
