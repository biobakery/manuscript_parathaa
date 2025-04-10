## Benchmarks full length assignments
require(docopt)

'Usage:
  full_length_bench.R [-p <parathaa_PATH> --dada_db_FL <dada2_db> -t <input_taxonomy> -o <output> --paraAssign <para_Assignments> --query <query seqs> -s <seed_data> -b <min boot>]
  
  
Options:
  -p directory where parathaa github repo is cloned
  --dada_db_FL location of DADA2 db
  --paraAssign location of parathaa assignments
  --query location of query reads 
  -t input taxonomy
  -o output
  -s seed_data
  -b minboot [default=80]
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


source("src/performance.table.R")


## This is very similar to the variable region function with some exceptions to deal with full length sequences
run.full.bench <- function(parathaaFile, sequenceFile, outputDir,
                               DADAdb, inFileTaxdata, inFileSeedDB, minboot=80, historic=FALSE){
  
  dir.create(outputDir, recursive = T, showWarnings = F)

  message(parathaaFile)
  # Read in parathaa data
  parathaaData <- read.delim(
    parathaaFile,
    sep='\t', fill=T, stringsAsFactors = F, header=T)
  message("DONE READING PARATHAA")
  #select only taxonomic data and group by the query.name
  tax_parathaa <- parathaaData %>%
    dplyr::select(query.name, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    group_by(query.name) 
  
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
  colnames(otumat) <- paste0("Parathaa.", "FL")
  OTU_parathaa <- otu_table(otumat, taxa_are_rows = TRUE)
  
  #create sample information so it can be imported into phyloseq
  samp_parathaa <- data.frame(colnames(OTU_parathaa), rep("parathaa", length(colnames(OTU_parathaa))),
                              rep("FL", length(colnames(OTU_parathaa))))
  colnames(samp_parathaa) <- c("sampleID", "Taxonomy_type", "Region")
  rownames(samp_parathaa) <- samp_parathaa$sampleID
  
  #print out data about the sample data frame
  str(samp_parathaa)
  SAMP_parathaa <- sample_data(samp_parathaa)
  
  #combined into phyloseq object and print out its basic information
  ps1_parathaa <- phyloseq(OTU_parathaa, TAX_parathaa, SAMP_parathaa)
  print(ps1_parathaa)
  
  
  
  
  ## Assign taxonomy with DADA2
  getNames <- read.fasta(file(sequenceFile), as.string = TRUE,
                         forceDNAtolower = FALSE, whole.header = FALSE)
  #get the names and split it by (tab as we only want to keep the accession ID)
  names1 <- str_split(getName(getNames), "\t", simplify=TRUE)
  #remove > in fasta headers if they exist
  names1 <- names1[,1] %>%
    str_remove(">")
  #create data frame with sequence and ID
  name.df <- data.frame("sequence" = unlist(getSequence(getNames, as.string=T)), taxaIDs = names1)
  nChars2 <- name.df %>% filter(str_detect(sequence, "N|M|R|K|Y|S|W|D|B|H|V")) %>% pull(taxaIDs)
  
  taxa_dada2 <- assignTaxonomy(sequenceFile, 
                         DADAdb,
                         multithread=TRUE,
                         minBoot = minboot)

  stopifnot(identical(tolower(name.df$sequence), tolower(unname(rownames(taxa_dada2)))))
  rownames(taxa_dada2) <- name.df$taxaIDs
  #select so we only keep taxonomy columns
  taxa_dada2 <- as.data.frame(taxa_dada2)
  taxa_dada2$taxaIDs <- rownames(taxa_dada2)
  
  tax_dada3 <- taxa_dada2 %>%
    select(Kingdom, Phylum, Class, Order, Family, Genus, Species) 
  
  #Add Species name so that genus is also there and deal with allowing multi assignments so its the same 
  # format as Parathaa:
  tax_dada3 <- as.matrix(tax_dada3)
  ## Place DADA2 taxonomies into phyloseq object
  TAX_dada <- tax_table(tax_dada3)
  
  ## Make dumby OTU table for the phyloseq object
  otutab <- as.data.frame(taxa_dada2) %>% 
    select(taxaIDs, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    dplyr::group_by(taxaIDs, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
    dplyr::count(taxaIDs, name="DADA2") 
  otumat <- as.matrix(otutab$DADA2)
  rownames(otumat) <- otutab$taxaIDs
  colnames(otumat) <- paste0("DADA2.", "FL")
  OTU_dada <- otu_table(otumat, taxa_are_rows = TRUE)
  
  ## Make sample data for the phyloseq object
  samp_dada <- data.frame(colnames(OTU_dada), rep("DADA2", length(colnames(OTU_dada))),
                          rep("FL", length(colnames(OTU_dada))))
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
  
  if(SILVA){
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
  }else{
    taxdata <- taxdata %>%
      tidyr::separate(path, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep="\\|", fill="right")
    
  }

  
  
  ## Read in seed db to identify what taxonomy the classifiers are actually aware of

  
  if(SILVA){
    SeedTax <- read.table(inFileSeedDB , header=F, fill=TRUE,sep='\t')
    SeedTax <- SeedTax %>%
      separate(col=V1, into=c("primaryAccession", "ArbID"), sep="\\.") %>%
      separate(col=V2, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep=";", extra="drop") %>%
      filter(Kingdom=="Bacteria" & !is.na(Genus) & Genus!="") 
  }else{
    SeedTax <- read.table(inFileSeedDB , header=T, fill=TRUE,sep='\t')
    SeedTax <- SeedTax %>% tidyr::separate(col=path, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Spec"), sep="\\|")
  }
  


  
  
  
  ## Make synthetic comparison dataset
  
  # Grab the parathaa taxonomy table
  synth.parathaa<- as.data.frame(tax_table(ps1_parathaa))

  
  #Join the Parathaa taxonomy table to the true reference taxonomies by their accession
  if(SILVA){
    synth.parathaa$AccID <- rownames(synth.parathaa)
    synth.parathaa2 <- left_join(synth.parathaa, taxdata, by="AccID")
  }else{
    synth.parathaa$primaryAccession <- rownames(synth.parathaa)
    synth.parathaa2 <- left_join(synth.parathaa, taxdata, by="primaryAccession")
  }
  
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
  if(SILVA){
    synth.dada$AccID <- rownames(synth.dada)
    synth.dada2 <- left_join(synth.dada, taxdata, by="AccID")
  }else{
    synth.dada$primaryAccession <- rownames(synth.dada)
    synth.dada2 <- left_join(synth.dada, taxdata, by="primaryAccession")
  }

  
  #Comapre the reference taxonomy to dada2 assigned taxonomy in the same manner as above

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
  
  colnames(synth.dada3) <- gsub("\\.x", ".dada", colnames(synth.dada3))
  colnames(synth.dada3) <- gsub("\\.y", ".silva", colnames(synth.dada3)) 
  
  colnames(synth.parathaa3) <- gsub("\\.x", ".parathaa", colnames(synth.parathaa3))
  colnames(synth.parathaa3) <- gsub("\\.y", ".silva", colnames(synth.parathaa3))
  
  
  #Join the two comparison dataframes
  if(SILVA){
    compare.synth <- dplyr::full_join(synth.dada3, synth.parathaa3, by="AccID")
  }else{
    compare.synth <- dplyr::full_join(synth.dada3, synth.parathaa3, by="primaryAccession")
    
  }
  #Rename columns to what tool they correspond to the tools, note if a different tool is used we still use this...
  if(!identical(compare.synth$Species.silva.x, compare.synth$Species.silva.y)){
    #fill in missing taxa in case some sequences were filtered out by either tool.
    missing_taxonomy_x <- which(is.na(compare.synth$Species.silva.x))
    
    #this basically protects against cases where the column we use for Species.silva is not totally filled out
    #as it should be.
    stopifnot(length(missing_taxonomy_x)==0)
  }
  
  compare.synth <- compare.synth %>% 
    mutate(Species.silva = word(Species.silva.x, 1, 2)) %>%
    rename(Genus.silva = Genus.silva.x)
  
  # Remove seqs with N characters: 
  if(length(nChars2) > 0)
    compare.synth <- compare.synth %>% filter(!AccID %in% nChars2)
  
  # save the comparison data
  save(compare.synth, file= file.path(outputDir, paste0("FL", "_full_comparisons.RData")))
  
  #write overall performance metrics for the comparison
  t1 <- performance.table(compare.synth, "Species")
  t1 <- t1[-which(rownames(t1)=="Unassigned Correct"),] 
  rownames(t1)[which(rownames(t1)=="Unassigned Incorrect")] <- "Unassigned"
  write.table(t1, file= file.path(outputDir, paste0("FL", "_Species_performance.tsv")), sep="\t", col.names = NA)
  t2 <- performance.table(compare.synth, "Genus")
  t2 <- t2[-which(rownames(t2)=="Unassigned Correct"),] 
  rownames(t2)[which(rownames(t2)=="Unassigned Incorrect")] <- "Unassigned"
  write.table(t2, file= file.path(outputDir, paste0("FL", "_Genus_performance.tsv")), sep="\t", col.names = NA)
  
  
  ## We next do the same as above but with an adjusted metric
  ## In the adjusted metric we treat unassignments as "correct" if the underlying seed database didn't have that
  ## taxonomic group within it
  
  #we only run this if using SILVA we don't need to do this for GTDB benching since we already filtered the data
  #so that sequences that are used are represented in the old GTDB database.
  if(!historic){
    ## would need to update this line since accessions are not consistent between GTDB releases 
    ## if we want to bench historically.
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
    write.table(t3, file= file.path(outputDir, paste0("FL", "_Species_performance_adjust.tsv")), sep="\t", col.names = NA)
    t4 <- performance.table(compare.synth_adjust, "Genus")
    write.table(t4, file= file.path(outputDir, paste0("FL", "_Genus_performance_adjust.tsv")), sep="\t", col.names = NA)
    
  }

}



run.full.bench(parathaaFile = opts$paraAssign, 
               sequenceFile = opts$query, 
               outputDir = opts$o, 
               DADAdb = opts$dada_db_FL, 
               inFileTaxdata = opts$t, 
               inFileSeedDB = opts$s,
               minboot = as.numeric(opts$b)
)


