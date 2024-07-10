## Compare taxonomies from DADA2 vs PARATHAA


require(docopt)

'Usage:
  Plots.Table.1.R [-p <parathaa_PATH> --dada_db <dada2_db> --dada_db_sp <dada2_db_sp> -t <input_taxonomy> -o <output> --paraAssignV4V5 <V4V5_taxonomy_file> --paraAssignV1V2 <V1V2_taxonomy_file> --queryV4V5 <V4V5 seqs> --queryV1V2 <V1V2 seqs> -s <seed_data>]
  
  
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
  -s seed_data
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

##There is probably a better way of sourcing this file... but for now leave it as is...
source("src/performance.table.R")



##################################
## Synthetic community analysis ##
##################################


## Define Function ##
run.synthetic.data <- function(parathaaFile, sequenceFile, regionName, outputDir, dadaAllowMult = FALSE,
                               DADAdb, DADAdb.sp, inFileTaxdata, inFileSeedDB){
  
  dir.create(outputDir, recursive = T, showWarnings = F)
  
  # Read in parathaa data
  parathaaData <- read.delim(
    parathaaFile,
    sep='\t', fill=T, stringsAsFactors = F, header=T)
  
  #select only taxonomic data and group by the query.name
  tax_parathaa <- parathaaData %>%
    dplyr::select(query.name, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    group_by(query.name) 
  
  ## we treat unclassified labels as unassigned in this case
  hierarchy <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  tax_parathaa <- tax_parathaa %>% mutate_at(vars(hierarchy), ~ str_replace(., "\\b\\w+ Unclassified", ""))
  ##then if we have blank string we replace with NA
  tax_parathaa <- tax_parathaa %>% mutate_at(vars(hierarchy), ~ na_if(., ""))
  ## replace any that are just ;
  tax_parathaa <- tax_parathaa %>% mutate_at(vars(hierarchy), ~ str_replace(., "^;+$", "")) 
  
  #convert to a matrix so that it can be input into a phyloseq tax_table
  taxmat <- tax_parathaa %>% as.matrix
  rownames(taxmat) <- tax_parathaa$query.name
  taxmat <- taxmat[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]  
  TAX_parathaa <- tax_table(taxmat)
  
  #create phyloseq otu table from parathaa data with a dumby count column named Parathaa
  otutab <- tax_parathaa %>% 
    dplyr::group_by(query.name, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
    dplyr::count(query.name, name="Parathaa")
  #convert to a matrix
  otumat <- as.matrix(otutab$Parathaa)
  rownames(otumat) <- otutab$query.name
  colnames(otumat) <- paste0("Parathaa.", regionName)
  OTU_parathaa <- otu_table(otumat, taxa_are_rows = TRUE)
  
  #create sample information so it can be imported into phyloseq
  samp_parathaa <- data.frame(colnames(OTU_parathaa), rep("parathaa", length(colnames(OTU_parathaa))),
                              rep(regionName, length(colnames(OTU_parathaa))))
  colnames(samp_parathaa) <- c("sampleID", "Taxonomy_type", "Region")
  rownames(samp_parathaa) <- samp_parathaa$sampleID
  
  #print out data about the sample data frame
  str(samp_parathaa)
  SAMP_parathaa <- sample_data(samp_parathaa)
  
  #combined into phyloseq object and print out its basic information
  ps1_parathaa <- phyloseq(OTU_parathaa, TAX_parathaa, SAMP_parathaa)
  print(ps1_parathaa)
  
  ## Assign taxonomy with DADA2
  
  ## First, get names and sequences from fasta file
  getNames <- read.fasta(file(sequenceFile), as.string = TRUE,
                         forceDNAtolower = FALSE, whole.header = FALSE)
  #get the names and split it by (tab as we only want to keep the accession ID)
  names1 <- str_split(getName(getNames), "\t", simplify=TRUE)
  #remove > in fasta headers if they exist
  names1 <- names1[,1] %>%
    str_remove(">")
  #create data frame with sequence and ID
  name.df <- data.frame("sequence" = unlist(getSequence(getNames, as.string=T)), taxaIDs = names1)
  
  #grab sequences with N in them as they need to be removed for species assignment by DADA2
  nChars2 <- name.df %>% filter(str_detect(sequence, "N|M|R|K|Y|S|W|D|B|H")) %>% pull(taxaIDs)
  
  ## Next, assign taxonomy to genus level with DADA2 (takes a few minutes)
  set.seed(3874)
  taxa <- assignTaxonomy(sequenceFile, 
                         DADAdb,
                         multithread=TRUE)
  
  
  ## Remove sequences with undefined ("N") bases, store until after species assignment
  taxa.test <- as.data.frame(taxa)
  taxa.test$taxaIDs <- names1
  nChars <- grep("N|M|R|K|Y|S|W|D|B|H", rownames(taxa.test))
  print(paste("Removing", length(nChars), "sequences with N bases"))
  withNbases <- taxa.test[nChars,]
  taxa <- taxa[-nChars,]
  
  ## Perform species assignment with DADA2 (takes a few minutes)
  taxa.sp <- addSpecies(taxa, DADAdb.sp, allowMultiple = dadaAllowMult)
  
  
  ## Add in reference IDs and taxonomy from sequences with "N" bases
  tax_dada <- as.data.frame(taxa.sp)
  tax_dada$sequence <- str_split(rownames(tax_dada), "\\.", simplify=TRUE)[,1]
  getnamSubset <-name.df %>% filter(sequence %in% tax_dada$sequence)
  tax_dada2 <- cbind(tax_dada, "taxaIDs" =getnamSubset$taxaIDs)
  tax_dada2 <- full_join(tax_dada2, withNbases)
  rownames(tax_dada2) <- tax_dada2[,"taxaIDs"]
  
  #select so we only keep taxonomy columns
  tax_dada3 <- tax_dada2 %>%
    select(Kingdom, Phylum, Class, Order, Family, Genus, Species) 
  
  #Add Species name so that genus is also there and deal with allowing multi assignments so its the same 
  # format as Parathaa:
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
  
  ## Make dumby OTU table for the phyloseq object
  otutab <- as.data.frame(tax_dada2) %>% 
    select(taxaIDs, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    dplyr::group_by(taxaIDs, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
    dplyr::count(taxaIDs, name="DADA2") 
  otumat <- as.matrix(otutab$DADA2)
  rownames(otumat) <- otutab$taxaIDs
  colnames(otumat) <- paste0("DADA2.", regionName)
  OTU_dada <- otu_table(otumat, taxa_are_rows = TRUE)
  
  ## Make sample data for the phyloseq object
  samp_dada <- data.frame(colnames(OTU_dada), rep("DADA2", length(colnames(OTU_dada))),
                          rep(regionName, length(colnames(OTU_dada))))
  colnames(samp_dada) <- c("sampleID", "Taxonomy_type", "Region")
  rownames(samp_dada) <- samp_dada$sampleID
  #print out sample data for phyloseq object
  str(samp_dada)
  SAMP_dada <- sample_data(samp_dada)
  
  ps1_dada <- phyloseq(OTU_dada, TAX_dada, SAMP_dada)
  #print dada2 phyloseq object
  print(ps1_dada)
  
  
  ############################################
  ### Assess performance on Synthetic Data ###
  ############################################
  
  #Get reference taxdata from SILVA:
  taxdata <- read.table(inFileTaxdata , header=T, fill=TRUE,sep='\t', quote="")
  taxdata <- taxdata %>%
    unite("AccID", c("primaryAccession", "start", "stop"), sep=".", remove=F)
  taxdata <- taxdata %>%
    mutate(taxonomy=paste0(path, organism_name))
  
  #remove Eukaryota Kingdom as they have taxonomy that we are not handling here.
  taxdata <- taxdata %>% filter(!grepl("^Eukaryota;", path))
  
  taxdata <- taxdata %>%
    select(AccID, primaryAccession, start, stop, taxonomy) %>%
    separate(col=taxonomy, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";", fill="right")
  
  ## Fix up remaining silva taxonomy by removing sub species and 'uncultured' species
  taxdata <- SILVA.species.editor(taxdata)
  
  
  ## Read in seed db to identify what taxonomy the classifiers are actually aware of
  SeedTax <- read.table(inFileSeedDB , header=F, fill=TRUE,sep='\t')
  
  SeedTax <- SeedTax %>%
    separate(col=V1, into=c("primaryAccession", "ArbID"), sep="\\.") %>%
    separate(col=V2, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep=";", extra="drop") %>%
    filter(Kingdom=="Bacteria" & !is.na(Genus) & Genus!="") 
  
  ### Look over this code to figure out why we are doing it this way ###
  taxdata_seedless <- taxdata %>% 
    filter(!primaryAccession %in% SeedTax$primaryAccession) %>%
    ## Subset to Genera in Seed for now
    ## Come back to why we need to do this?
    filter(Genus %in% unique(SeedTax$Genus))
  

  
  ## Make synthetic comparison dataset
  
  # Grab the parathaa taxonomy table
  synth.parathaa<- as.data.frame(tax_table(ps1_parathaa))
  synth.parathaa$AccID <- rownames(synth.parathaa)
  
  #Join the Parathaa taxonomy table to the true reference taxonomies by their accession
  synth.parathaa2 <- left_join(synth.parathaa, taxdata, by="AccID")
  
  #Format species name to match references
  synth.parathaa2 <- synth.parathaa2 %>% 
    mutate(Species.x = unlist(lapply(str_split(Species.x, ";"), FUN=function(x) paste0(word(x,1,2), collapse = ";" ))))
  #set NAs as either NA strings or empty strings
  synth.parathaa2 <- synth.parathaa2 %>% 
    mutate(Species.x = ifelse(Species.x=="NA", NA, Species.x),
           Genus.x = ifelse(Genus.x=="" | Genus.x=="uncultured", NA, Genus.x))
  #Compare the reference taxonomy to assigned taxonomy by parathaa
  #If species are the same set Flag to true
  #If genera are the same set Flag.genus to true
  #In both cases we call a true match if at least one matchs when there are multiple assignments
  synth.parathaa3 <- synth.parathaa2 %>% 
    dplyr::rowwise() %>%
    mutate(Flag = ifelse(is.na(Species.x), NA, word(Species.y, 1, 2) %in% str_split(Species.x, ";", simplify = T)),
           Flag.genus = ifelse(is.na(Genus.x), NA, Genus.y %in% str_split(Genus.x, ";", simplify = T))
    )
  
  #Grab dada2 taxonomy table
  synth.dada <- as.data.frame(tax_table(ps1_dada))
  synth.dada$AccID <- rownames(synth.dada)
  
  #Comapre the reference taxonomy to dada2 assigned taxonomy in the same manner as above
  synth.dada2 <- left_join(synth.dada, taxdata, by="AccID")
  synth.dada2 <- synth.dada2 %>% 
    mutate(Species.x = unlist(lapply(str_split(Species.x, ";"), FUN=function(x) paste0(word(x,1,2), collapse = ";" ))))
  synth.dada2 <- synth.dada2 %>% 
    mutate(Species.x = ifelse(Species.x=="NA", NA, Species.x),
           Genus.x = ifelse(Genus.x=="" | Genus.x=="uncultured", NA, Genus.x))
  synth.dada3 <- synth.dada2 %>% 
    dplyr::rowwise() %>%
    mutate(Flag = ifelse(is.na(Species.x), NA, word(Species.y, 1, 2) %in% str_split(Species.x, ";", simplify = T)),
           Flag.genus = ifelse(is.na(Genus.x), NA, Genus.y %in% str_split(Genus.x, ";", simplify = T))
    )
  
  
  #Join the two comparison dataframes
  compare.synth <- dplyr::full_join(synth.dada3, synth.parathaa3, by="AccID")
  #Rename columns to what tool they correspond to
  compare.synth <- compare.synth %>% 
    mutate(Species.silva = word(Species.y.y, 1, 2)) %>%
    rename(Species.parathaa = Species.x.y,
           Species.dada = Species.x.x,
           Genus.dada = Genus.x.x,
           Genus.parathaa = Genus.x.y,
           Genus.silva = Genus.y.y)
  
  # Remove seqs with N characters: 
  compare.synth <- compare.synth %>% filter(!AccID %in% nChars2)
  # save the comparison data
  save(compare.synth, file= file.path(outputDir, paste0(regionName, "_full_comparisons.RData")))
  
  #write overall performance metrics for the comparison
  t1 <- performance.table(compare.synth, "Species")
  write.table(t1, file= file.path(outputDir, paste0(regionName, "_Species_performance.tsv")), sep="\t", col.names = NA)
  t2 <- performance.table(compare.synth, "Genus")
  write.table(t2, file= file.path(outputDir, paste0(regionName, "_Genus_performance.tsv")), sep="\t", col.names = NA)
  
  
  ## We next do the same as above but with an adjusted metric
  ## In the adjusted metric we treat unassignments as "correct" if the underlying seed database didn't have that
  ## taxonomic group within it
  
  ###Seed_genus is TRUE if the genus was included in the seed database. We do not want to change these assignments..
  taxdata_seed <- taxdata %>% filter(primaryAccession %in% SeedTax$primaryAccession)
  taxdata_SP <- word(taxdata_seed$Species, 1, 2)
  taxdata_SP <- taxdata_SP[-which(is.na(taxdata_SP))]
  compare.synth_adjust <- compare.synth %>% mutate(seed_genus=Genus.silva %in% taxdata_seed$Genus)
  compare.synth_adjust <- compare.synth_adjust %>% mutate(seed_species=Species.silva %in% taxdata_SP)
  
  #Make a corrected genus flag where we set previous NAs to true if they are unassigned and not in the seed DB
  compare.synth_adjust <- compare.synth_adjust %>% mutate(Flag.genus.x_cor=ifelse(is.na(Flag.genus.x) & !seed_genus, TRUE, Flag.genus.x))
  compare.synth_adjust <- compare.synth_adjust %>% mutate(Flag.genus.y_cor=ifelse(is.na(Flag.genus.y) & !seed_genus, TRUE, Flag.genus.y))
  
  #Same as above for species
  compare.synth_adjust <- compare.synth_adjust %>% mutate(Flag.x_cor=ifelse(is.na(Flag.x) & !seed_species, TRUE, Flag.x))
  compare.synth_adjust <- compare.synth_adjust %>% mutate(Flag.y_cor=ifelse(is.na(Flag.y) & !seed_species, TRUE, Flag.y))
  
  
  compare.synth_adjust$Flag.genus.x <- compare.synth_adjust$Flag.genus.x_cor
  compare.synth_adjust$Flag.genus.y <- compare.synth_adjust$Flag.genus.y_cor
  
  compare.synth_adjust$Flag.x <- compare.synth_adjust$Flag.x_cor
  compare.synth_adjust$Flag.y <- compare.synth_adjust$Flag.y_cor
  
  t3 <- performance.table(compare.synth_adjust, "Species")
  write.table(t3, file= file.path(outputDir, paste0(regionName, "_Species_performance_adjust.tsv")), sep="\t", col.names = NA)
  t4 <- performance.table(compare.synth_adjust, "Genus")
  write.table(t4, file= file.path(outputDir, paste0(regionName, "_Genus_performance_adjust.tsv")), sep="\t", col.names = NA)
  
  save(compare.synth_adjust, file= file.path(outputDir, paste0(regionName, "_full_comparisons_adjust.RData")))
}

## Define variables used across all function calls
## 
DADAdb <- opts$dada_db
DADAdb.sp <- opts$dada_db_sp
inFileTaxdata <- opts$t
inFileSeedDB <- opts$s


##### CALL FUNCTION ######
run.synthetic.data(parathaaFile = opts$paraAssignV4V5, 
                   sequenceFile = opts$queryV4V5,
                   regionName = "V4V5", 
                   outputDir=paste(opts$o, "/Figures/synth_mult_arc/", sep=""),
                   dadaAllowMult = T,
                   DADAdb = DADAdb,
                   DADAdb.sp = DADAdb.sp,
                   inFileTaxdata = inFileTaxdata,
                   inFileSeedDB = inFileSeedDB)

run.synthetic.data(parathaaFile = opts$paraAssignV1V2, 
                   sequenceFile = opts$queryV1V2, 
                   regionName = "V1V2", 
                   outputDir=paste(opts$o, "/Figures/synth_mult_arc", sep=""),
                   dadaAllowMult = T,
                   DADAdb = DADAdb,
                   DADAdb.sp = DADAdb.sp,
                   inFileTaxdata = inFileTaxdata,
                   inFileSeedDB = inFileSeedDB)

run.synthetic.data(parathaaFile = opts$paraAssignV4V5, 
                   sequenceFile = opts$queryV4V5,
                   regionName = "V4V5", 
                   outputDir=paste(opts$o, "/Figures/synth_nomult_arc/", sep=""),
                   dadaAllowMult = F,
                   DADAdb = DADAdb,
                   DADAdb.sp = DADAdb.sp,
                   inFileTaxdata = inFileTaxdata,
                   inFileSeedDB = inFileSeedDB)

run.synthetic.data(parathaaFile = opts$paraAssignV1V2, 
                   sequenceFile = opts$queryV1V2, 
                   regionName = "V1V2", 
                   outputDir=paste(opts$o, "/Figures/synth_nomult_arc", sep=""),
                   dadaAllowMult = F,
                   DADAdb = DADAdb,
                   DADAdb.sp = DADAdb.sp,
                   inFileTaxdata = inFileTaxdata,
                   inFileSeedDB = inFileSeedDB)

