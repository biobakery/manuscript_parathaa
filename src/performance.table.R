
## Function that runs the benchmarks against DADA2
run.synthetic.data <- function(parathaaFile, sequenceFile, regionName, outputDir, dadaAllowMult = FALSE,
                               DADAdb, DADAdb.sp, inFileTaxdata, inFileSeedDB, SILVA=TRUE, historic=FALSE,
                               full_length=FALSE, minboot=80){
  
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
  tax_parathaa <- tax_parathaa %>% mutate_at(vars(hierarchy), ~ str_replace(., "\\b\\w+ unclassified", ""))
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
  nChars2 <- name.df %>% filter(str_detect(sequence, "N|M|R|K|Y|S|W|D|B|H|V")) %>% pull(taxaIDs)
  
  ## Next, assign taxonomy to genus level with DADA2 (takes a few minutes)
  if(full_length){
    message(paste0("assigning taxonomy with minboot= ", minboot))
    tax_dada2 <- assignTaxonomy(sequenceFile, 
                                 DADAdb,
                                 multithread=TRUE,
                                 minBoot = minboot, 
                                outputBootstraps=T)
    save(tax_dada2, file= file.path(outputDir, paste0(regionName, "DADA2_assignments_with_confidence.RData")))
    tax_dada2 <- tax_dada2[[1]]
    
    stopifnot(identical(tolower(name.df$sequence), tolower(unname(rownames(tax_dada2)))))
    rownames(tax_dada2) <- name.df$taxaIDs
    #select so we only keep taxonomy columns
    tax_dada2 <- as.data.frame(tax_dada2)
    tax_dada2$taxaIDs <- rownames(tax_dada2)
    
    #benching against full length + exact matching
    
    
    #select so we only keep taxonomy columns
    tax_dada3 <- tax_dada2 %>%
      select(Kingdom, Phylum, Class, Order, Family, Genus, Species) 
    
    #Add Species name so that genus is also there and deal with allowing multi assignments so its the same 
    # format as Parathaa:
    tax_dada3 <- tax_dada3 %>% 
      as.matrix()
    
    rownames(tax_dada3) <- tax_dada2[,"taxaIDs"]
  }else{
    set.seed(3874)
    taxa <- assignTaxonomy(sequenceFile, 
                           DADAdb,
                           multithread=TRUE,
                           outputBootstraps = T)
    
    save(taxa, file= file.path(outputDir, paste0(regionName, "DADA2_assignments_with_confidence.RData")))
    
    taxa <- taxa[[1]]
    ## Remove sequences with undefined ("N") bases, store until after species assignment
    taxa.test <- as.data.frame(taxa)
    taxa.test$taxaIDs <- names1
    nChars <- grep("N|M|R|K|Y|S|W|D|B|H|V", rownames(taxa.test))
    print(paste("Removing", length(nChars), "sequences with N bases"))
    withNbases <- taxa.test[nChars,]
    if(length(nChars!=0))
      taxa <- taxa[-nChars,]
    
    ## Perform species assignment with DADA2 (takes a few minutes)
    taxa.sp <- addSpecies(taxa, DADAdb.sp, allowMultiple = dadaAllowMult)
    
    
    ## Add in reference IDs and taxonomy from sequences with "N" bases
    tax_dada <- as.data.frame(taxa.sp)
    tax_dada$sequence <- str_split(rownames(tax_dada), "\\.", simplify=TRUE)[,1]
    getnamSubset <-name.df %>% filter(tolower(sequence) %in% tolower(tax_dada$sequence))
    tax_dada2 <- cbind(tax_dada, "taxaIDs" =getnamSubset$taxaIDs)
    if(length(nChars!=0))
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
  }
  
  

  
  
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
    if(!historic){
      SeedTax <- read.table(inFileSeedDB , header=T, fill=TRUE,sep='\t')
      SeedTax <- SeedTax %>% tidyr::separate(col=path, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Spec"), sep="\\|")
    }
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
  
  ##rename columns so they don't get confused.
  
  colnames(synth.dada3) <- gsub("\\.x", ".dada", colnames(synth.dada3))
  colnames(synth.dada3) <- gsub("\\.y", ".silva", colnames(synth.dada3)) 
  
  colnames(synth.parathaa3) <- gsub("\\.x", ".parathaa", colnames(synth.parathaa3))
  colnames(synth.parathaa3) <- gsub("\\.y", ".silva", colnames(synth.parathaa3))
  
  #Join the two comparison dataframes
  if(SILVA){
    colnames(name.df)[2] <- "AccID"
    
    #dada info.
    compare.synth <- name.df %>% 
      left_join(synth.dada3[, grepl("dada|Flag|AccID", colnames(synth.dada3))], by="AccID")
    
    #add parathaa info
    compare.synth <- compare.synth %>% 
      left_join(synth.parathaa3[, grepl("parathaa|Flag|AccID", colnames(synth.parathaa3))], by="AccID")
    
    compare.synth <- compare.synth %>% 
      left_join(taxdata, by="AccID")
  }else{
    colnames(name.df)[2] <- "primaryAccession"
    
    #add dada assignments to all sequences that were attempted to assign profiles
    compare.synth <- name.df %>% 
      left_join(synth.dada3[, grepl("dada|Flag|primary", colnames(synth.dada3))], by="primaryAccession")
    #add parathaa assignments
    compare.synth <- compare.synth %>% 
      left_join(synth.parathaa3[, grepl("parathaa|Flag|primary", colnames(synth.parathaa3))], by="primaryAccession")
    compare.synth <- compare.synth %>% 
      left_join(taxdata, by="primaryAccession")
  }
  
  compare.synth <- compare.synth %>% 
    #some data in silva has no speices name
    mutate(Species.groundtruth = word(Species, 1, 2)) %>%
    rename(Genus.groundtruth = Genus) %>%
    rename(Family.groundtruth = Family) %>%
    rename(Order.groundtruth = Order) %>%
    rename(Class.groundtruth = Class) %>%
    rename(Phylum.groundtruth = Phylum) %>%
    rename(Kingdom.groundtruth = Kingdom) %>%
    select(-Species)
  
  
  # Remove seqs with N characters: 
  if(length(nChars2) > 0)
    compare.synth <- compare.synth %>% filter(!AccID %in% nChars2)
  # save the comparison data
  save(compare.synth, file= file.path(outputDir, paste0(regionName, "_full_comparisons.RData")))
  
  #write overall performance metrics for the comparison
  ##for species comparison we need to filter out taxa that don't have a ground truth speices
  t1 <- performance.table(compare.synth %>% filter(!is.na(Species.groundtruth)), "Species")
  t1 <- t1[-which(rownames(t1)=="Unassigned Correct"),] 
  rownames(t1)[which(rownames(t1)=="Unassigned Incorrect")] <- "Unassigned"
  write.table(t1, file= file.path(outputDir, paste0(regionName, "_Species_performance.tsv")), sep="\t", col.names = NA)
  
  ## Genus level results
  t2 <- performance.table(compare.synth %>% filter(!is.na(Genus.groundtruth)), "Genus")
  
  t2 <- t2[-which(rownames(t2)=="Unassigned Correct"),] 
  rownames(t2)[which(rownames(t2)=="Unassigned Incorrect")] <- "Unassigned"
  
  write.table(t2, file= file.path(outputDir, paste0(regionName, "_Genus_performance.tsv")), sep="\t", col.names = NA)
  
  
  ## We next do the same as above but with an adjusted metric
  ## In the adjusted metric we treat unassignments as "correct" if the underlying seed database didn't have that
  ## taxonomic group within it
  
  # We only need to do this for SILVA we don't need to do this for GTDB benching since we already filtered
  # the sequences we may want to do this in future if we do holdout GTDB benching but for now lets just
  # make this run only if SILVA is TRUE
  if(!historic){
    ## would need to update this line since accessions are not consistent between GTDB releases 
    ## if we want to bench historically.
    taxdata_seed <- taxdata %>% filter(primaryAccession %in% SeedTax$primaryAccession)
    taxdata_SP <- word(taxdata_seed$Species, 1, 2)
    
    if(SILVA){
      taxdata_SP <- taxdata_SP[-which(is.na(taxdata_SP))]
    }
    
    compare.synth_adjust <- compare.synth %>% mutate(seed_genus=Genus.groundtruth %in% taxdata_seed$Genus)
    compare.synth_adjust <- compare.synth_adjust %>% mutate(seed_species=Species.groundtruth %in% taxdata_SP)
    
    #Make a corrected genus flag where we set previous set to NA. Set to true if they are unassigned and not in the seed DB
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
  
}


## Function that computes the performance of the tool
performance.table <- function(compareData, level){
  
  #If we are comparing species level data
  if(level=="Species"){
    
    #We treat true positives as anything that we set the FLAG variable as true
    TP.parathaa <- compareData %>% filter(Flag.y) %>% nrow()
    #We consider false positives to be anything that the FLAG variable is set to false
    FP.parathaa <- compareData %>% filter(!Flag.y) %>% nrow()
    #We do not consider true negatives in our analysis
    TN.parathaa <- 0
    #False negatives is cases where taxonomy is not assigned and the flag variable is not set
    FN.parathaa <- compareData %>% filter(is.na(Species.parathaa) & is.na(Flag.y)) %>% nrow()
    
    #multi correct are cases where a species is set but it doesn't equal the silva species and the flag is true
    multCorrect.parathaa <- compareData %>% filter(!Species.parathaa==Species.groundtruth & Flag.y) %>% nrow() /nrow(compareData)
    #unique correct are cases when taxonomy match up
    uniqueCorrect.parathaa <- compareData %>% filter(Species.parathaa==Species.groundtruth) %>% nrow() /nrow(compareData)
    #unassigned correct are cases where we think not assigning a taxonomy is the correct choice
    #this occurs when the reference query is from a taxonomy that is outside of the database's knowledge
    unassignedCorrect.parathaa <- compareData %>% filter(is.na(Species.parathaa) & Flag.y) %>% nrow() /nrow(compareData)
    #Unassigned incorrect are cases when we don't make an assignment but the taxonomy we were suppose to assign does exist 
    # within the database
    unassignedIncorrect.parathaa <- compareData %>% filter(is.na(Species.parathaa) & is.na(Flag.y)) %>% nrow() /nrow(compareData)
    multIncorrect.parathaa <- compareData %>% filter(grepl(";", Species.parathaa)) %>% filter(!Flag.y) %>% nrow() /nrow(compareData)
    uniquelyIncorrect.parathaa <- compareData %>% filter(!grepl(";", Species.parathaa)) %>% filter(!Flag.y) %>% nrow() / nrow(compareData)
  }
  
  # For information about individual variables see above section
  if(level=="Genus"){
    TP.parathaa <- compareData %>% filter( Flag.genus.y) %>% nrow()
    FP.parathaa <- compareData %>% filter(!Flag.genus.y) %>% nrow()
    TN.parathaa <- 0
    FN.parathaa <- compareData %>% filter(is.na(Genus.parathaa) & is.na(Flag.genus.y)) %>% nrow()
    
    multCorrect.parathaa <- compareData %>% filter(!Genus.parathaa==Genus.groundtruth & Flag.genus.y) %>% nrow() /nrow(compareData)
    uniqueCorrect.parathaa <- compareData %>% filter(Genus.parathaa==Genus.groundtruth) %>% nrow() /nrow(compareData)
    unassignedCorrect.parathaa <- compareData %>% filter(is.na(Genus.parathaa) & Flag.genus.y) %>% nrow() /nrow(compareData)
    unassignedIncorrect.parathaa <- compareData %>% filter(is.na(Genus.parathaa) & is.na(Flag.genus.y)) %>% nrow() /nrow(compareData)
    multIncorrect.parathaa <- compareData %>% filter(grepl(";", Genus.parathaa)) %>% filter(!Flag.genus.y) %>% nrow() /nrow(compareData)
    uniquelyIncorrect.parathaa <- compareData %>% filter(!grepl(";", Genus.parathaa)) %>% filter(!Flag.genus.y) %>% nrow() / nrow(compareData)
  }
  
  #Calculate metrics
  accuracy.parathaa <- (TP.parathaa + TN.parathaa) / (TP.parathaa + TN.parathaa + FP.parathaa + FN.parathaa)
  precision.parathaa <- TP.parathaa / (TP.parathaa + FP.parathaa)
  recall.parathaa <- TP.parathaa / (TP.parathaa + FN.parathaa)
  f1.parathaa <- 2 * (precision.parathaa * recall.parathaa) / (precision.parathaa + recall.parathaa)
  fpr.parathaa <- FP.parathaa / nrow(compareData)
  
  #Do the same as above but for dada2
  if(level=="Species"){
    TP.dada <- compareData %>% filter( Flag.x) %>% nrow()
    FP.dada <- compareData %>% filter(!Flag.x) %>% nrow()
    TN.dada <- 0
    FN.dada <- compareData %>% filter(is.na(Species.dada) & is.na(Flag.x)) %>% nrow()
    multCorrect.dada <- compareData %>% filter(!Species.dada==Species.groundtruth & Flag.x) %>% nrow() /nrow(compareData)
    uniqueCorrect.dada <- compareData %>% filter(Species.dada==Species.groundtruth) %>% nrow() /nrow(compareData)
    unassignedCorrect.dada <- compareData %>% filter(is.na(Species.dada) & Flag.x) %>% nrow() /nrow(compareData)
    unassignedIncorrect.dada <- compareData %>% filter(is.na(Species.dada) & is.na(Flag.x)) %>% nrow() /nrow(compareData)
    multIncorrect.dada <- compareData %>% filter(grepl(";", Species.dada)) %>% filter(!Flag.x) %>% nrow() /nrow(compareData)
    uniquelyIncorrect.dada <- compareData %>% filter(!grepl(";", Species.dada)) %>% filter(!Flag.x) %>% nrow() / nrow(compareData)
  }
  if(level=="Genus"){
    TP.dada <- compareData %>% filter( Flag.genus.x) %>% nrow()
    FP.dada <- compareData %>% filter(!Flag.genus.x) %>% nrow()
    TN.dada <- 0
    FN.dada <- compareData %>% filter(is.na(Genus.dada) & is.na(Flag.genus.x)) %>% nrow()
    multCorrect.dada <- compareData %>% filter(!Genus.dada==Genus.groundtruth & Flag.genus.x) %>% nrow() /nrow(compareData)
    uniqueCorrect.dada <- compareData %>% filter(Genus.dada==Genus.groundtruth) %>% nrow() /nrow(compareData)
    unassignedCorrect.dada <- compareData %>% filter(is.na(Genus.dada) & Flag.genus.x) %>% nrow() /nrow(compareData)
    unassignedIncorrect.dada <- compareData %>% filter(is.na(Genus.dada) & is.na(Flag.genus.x)) %>% nrow() /nrow(compareData)
    multIncorrect.dada <- compareData %>% filter(grepl(";", Genus.dada)) %>% filter(!Flag.genus.x) %>% nrow() /nrow(compareData)
    uniquelyIncorrect.dada <- compareData %>% filter(!grepl(";", Genus.dada)) %>% filter(!Flag.genus.x) %>% nrow() / nrow(compareData)
  }
  
  accuracy.dada <- (TP.dada + TN.dada) / (TP.dada + TN.dada + FP.dada + FN.dada)
  precision.dada <- TP.dada / (TP.dada + FP.dada)
  recall.dada <- TP.dada / (TP.dada + FN.dada)
  f1.dada <- 2 * (precision.dada * recall.dada) / (precision.dada + recall.dada)
  fpr.dada <- FP.dada / nrow(compareData)

  
  #set up the output table
  rows1 <- c(             "Accuracy", "Precision", "Recall", "F1 Score",
                          "Uniquely Correct", "One-to-many Correct", "Incorrect", "Unassigned Correct", "Unassigned Incorrect", 
                          "One-to-many Incorrect", "Uniquely Incorrect")
  
  table.out <- matrix(NA, nrow = length(rows1), ncol=2)
  colnames(table.out) <- c("Parathaa", "DADA2")
  rownames(table.out) <- rows1
  table.out["Uniquely Correct", "Parathaa"] <- uniqueCorrect.parathaa
  table.out["Uniquely Correct", "DADA2"] <- uniqueCorrect.dada
  table.out["One-to-many Correct", "Parathaa"] <- multCorrect.parathaa
  table.out["One-to-many Correct", "DADA2"] <- multCorrect.dada
  table.out["Incorrect", "Parathaa"] <- fpr.parathaa
  table.out["Incorrect", "DADA2"] <- fpr.dada
  table.out["Unassigned Correct", "Parathaa"] <- unassignedCorrect.parathaa
  table.out["Unassigned Correct", "DADA2"] <- unassignedCorrect.dada
  table.out["Unassigned Incorrect", "Parathaa"] <- unassignedIncorrect.parathaa
  table.out["Unassigned Incorrect", "DADA2"] <- unassignedIncorrect.dada
  table.out["One-to-many Incorrect", "Parathaa"] <- multIncorrect.parathaa
  table.out["One-to-many Incorrect", "DADA2"] <- multIncorrect.dada
  table.out["Uniquely Incorrect", "Parathaa"] <- uniquelyIncorrect.parathaa
  table.out["Uniquely Incorrect", "DADA2"] <- uniquelyIncorrect.dada
  
  table.out["Accuracy", "Parathaa"] <- accuracy.parathaa
  table.out["Accuracy", "DADA2"] <- accuracy.dada
  table.out["Precision", "Parathaa"] <- precision.parathaa
  table.out["Precision", "DADA2"] <- precision.dada
  table.out["Recall", "Parathaa"] <- recall.parathaa
  table.out["Recall", "DADA2"] <- recall.dada
  table.out["F1 Score", "Parathaa"] <- f1.parathaa
  table.out["F1 Score", "DADA2"] <- f1.dada


  print(round(table.out, 3))
}



Bench_IDTAXA <- function(sequenceFile, genus_class, species_class, threads=8, inFileTaxdata, 
                         inFileSeedDB, SILVA=T, outputDir, historic=FALSE, regionName){
  library(dada2)
  library(seqinr)
  
  ## load 

  
  getNames <- seqinr::read.fasta(file(sequenceFile), as.string = TRUE,
                         forceDNAtolower = FALSE, whole.header = FALSE)
  #get the names and split it by (tab as we only want to keep the accession ID)
  names1 <- str_split(getName(getNames), "\t", simplify=TRUE)
  #remove > in fasta headers if they exist
  names1 <- names1[,1] %>%
    str_remove(">")
  #create data frame with sequence and ID
  name.df <- data.frame("sequence" = unlist(getSequence(getNames, as.string=T)), taxaIDs = names1)
  
  #grab sequences with N in them as they need to be removed for species assignment by DADA2
  nChars2 <- name.df %>% filter(str_detect(sequence, "N|M|R|K|Y|S|W|D|B|H|V")) %>% pull(taxaIDs)
  
  dna <- DNAStringSet(getSequences(sequenceFile))
  ##genus level assignments
  genus_classifer <- readRDS(genus_class)
  
  ids_genus <- IdTaxa(dna, genus_classifer, type = "extended", strand = "top", processors = threads, verbose = TRUE)  # the strands SHOULD match, but I guess you can run strand = "both" to be certain
  
  ##turn into table to be used
  genus_assignments <- data.frame(do.call(rbind, ids_genus))
  if(regionName=="FL"){
    genus_assignments$AccID <- gsub(" .*", "", rownames(genus_assignments))
  }else{
    genus_assignments$AccID <- gsub("\t.*", "", rownames(genus_assignments))
  }

  
  genus_assignments <- genus_assignments %>% mutate(confidence_level=map(genus_assignments$confidence, function(x) x[1:7])) %>%
    unnest_wider(confidence_level, names_sep = "_")
  
  genus_assignments <- genus_assignments %>% unnest_wider(taxon, names_sep = "_")

  genus_assignments <- genus_assignments %>%
    dplyr::rename(Root_Assignment=taxon_1) %>%
    dplyr::rename(Kingdom_Assignment=taxon_2) %>%
    dplyr::rename(Phylum_Assignment=taxon_3) %>%
    dplyr::rename(Class_Assignment=taxon_4) %>%
    dplyr::rename(Order_Assignment=taxon_5) %>%
    dplyr::rename(Family_Assignment=taxon_6) %>%
    dplyr::rename(Genus_Assignment=taxon_7)
  
  genus_assignments$Kingdom_Assignment[grep("unclassified_", genus_assignments$Kingdom_Assignment, ignore.case = T)] <- NA
  genus_assignments$Phylum_Assignment[grep("unclassified_", genus_assignments$Phylum_Assignment, ignore.case = T)] <- NA
  genus_assignments$Class_Assignment[grep("unclassified_", genus_assignments$Class_Assignment, ignore.case = T)] <- NA
  genus_assignments$Order_Assignment[grep("unclassified_", genus_assignments$Order_Assignment, ignore.case = T)] <- NA
  genus_assignments$Family_Assignment[grep("unclassified_", genus_assignments$Family_Assignment, ignore.case = T)] <- NA
  genus_assignments$Genus_Assignment[grep("unclassified_", genus_assignments$Genus_Assignment, ignore.case = T)] <- NA

  ## create comparison frame
  if(regionName=="FL"){
    ids <- gsub(" .*", "", names(dna))
  }else{
    ids <- gsub("\t.*", "", names(dna))
  }

  
  #grab ground truths
  
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
    
    
    compare_frame <- data.frame("AccID"=ids)
    
    compare_frame <- compare_frame %>% left_join(genus_assignments)
    
    compare_frame <- compare_frame %>% left_join(taxdata)
    
    #set NAs as either NA strings or empty strings
    compare_frame2 <- compare_frame %>% 
      mutate(Genus = ifelse(Genus=="" | Genus=="uncultured", NA, Genus),
             Genus_Assignment = ifelse(Genus_Assignment=="" | Genus_Assignment=="uncultured", NA, Genus_Assignment))
    
    
    #Compare the reference taxonomy to assigned taxonomy by parathaa
    #If species are the same set Flag to true
    #If genera are the same set Flag.genus to true
    #In both cases we call a true match if at least one matchs when there are multiple assignments
    
    compare_frame3 <- compare_frame2 %>% 
      dplyr::rowwise() %>%
      #first we check if the assignment is Na if so we set genus flag to NA
      #next we check if the true genus is the assigned genus
      mutate(Flag.genus = ifelse(is.na(Genus_Assignment), NA, Genus %in% str_split(Genus_Assignment, ";", simplify = T)))
    
    # calculate performance for a single run and report
    ##remove Nchars
    if(length(nChars2) > 0)
      compare_frame3 <- compare_frame3 %>% filter(!AccID %in% nChars2)
    
    #no adjustment
    genus_performance <- performance_single_run(compare_frame3, level="Genus")
    genus_performance <- genus_performance[-which(rownames(genus_performance)=="Unassigned Correct"),,drop=FALSE] 
    #fix it
    rownames(genus_performance)[which(rownames(genus_performance)=="Unassigned Incorrect")] <- "Unassigned"
    
    write.table(genus_performance, file= file.path(outputDir, paste0(regionName, "_Genus_performance.tsv")), sep="\t", col.names = NA)
  }
  
  
  ## now get species level assignments
  specs_classifier <- readRDS(species_class)
  ids_specs <- IdTaxa(dna, specs_classifier, type = "extended", strand = "top", processors = threads, verbose = TRUE)  # the strands SHOULD match, but I guess you can run strand = "both" to be certain
  
  ##turn into table to be used
  specs_assignments <- data.frame(do.call(rbind, ids_specs))
  if(regionName=="FL"){
    specs_assignments$AccID <- gsub(" .*", "", rownames(specs_assignments))
  }else{
    specs_assignments$AccID <- gsub("\t.*", "", rownames(specs_assignments))
    
  }
  
  specs_assignments <- specs_assignments %>% mutate(confidence_level=map(specs_assignments$confidence, function(x) x[1:7])) %>%
    unnest_wider(confidence_level, names_sep = "_")
  
  specs_assignments <- specs_assignments %>% unnest_wider(taxon, names_sep = "_")
  
  specs_assignments <- specs_assignments %>%
    dplyr::rename(Root_Assignment=taxon_1) %>%
    dplyr::rename(Kingdom_Assignment=taxon_2) %>%
    dplyr::rename(Phylum_Assignment=taxon_3) %>%
    dplyr::rename(Class_Assignment=taxon_4) %>%
    dplyr::rename(Order_Assignment=taxon_5) %>%
    dplyr::rename(Family_Assignment=taxon_6) %>%
    dplyr::rename(Genus_Assignment=taxon_7) %>%
    dplyr::rename(Species_Assignment=taxon_8)
  
  #search through and set Na to all _unclassified
  specs_assignments$Kingdom_Assignment[grep("unclassified_", specs_assignments$Kingdom_Assignment, ignore.case = T)] <- NA
  specs_assignments$Phylum_Assignment[grep("unclassified_", specs_assignments$Phylum_Assignment, ignore.case = T)] <- NA
  specs_assignments$Class_Assignment[grep("unclassified_", specs_assignments$Class_Assignment, ignore.case = T)] <- NA
  specs_assignments$Order_Assignment[grep("unclassified_", specs_assignments$Order_Assignment, ignore.case = T)] <- NA
  specs_assignments$Family_Assignment[grep("unclassified_", specs_assignments$Family_Assignment, ignore.case = T)] <- NA
  specs_assignments$Genus_Assignment[grep("unclassified_", specs_assignments$Genus_Assignment, ignore.case = T)] <- NA
  specs_assignments$Species_Assignment[grep("unclassified_", specs_assignments$Species_Assignment, ignore.case = T)] <- NA
  
  specs_assignments <- specs_assignments %>% 
    mutate(Species_Assignment = word(Species_Assignment, 1, 2))
  
  
  if(SILVA){
    compare_frame_spec <- data.frame("AccID"=ids)
    
    compare_frame_spec <- compare_frame_spec %>% left_join(specs_assignments)
    
    compare_frame_spec <- compare_frame_spec %>% left_join(taxdata)
    
    #set ground truth to just genus species ignore any potential straing/sub-species etc.
    compare_frame_spec <- compare_frame_spec %>% 
      mutate(Species = word(Species, 1, 2))
    
    #set NAs as either NA strings or empty strings
    compare_frame3_spec <- compare_frame_spec %>% 
      rowwise() %>%
      mutate(Flag=ifelse(is.na(Species_Assignment), NA, word(Species, 1, 2) %in% str_split(Species_Assignment, ";", simplify = T)))

    # calculate performance for a single run and report
    ##remove Nchars
    if(length(nChars2) > 0)
      compare_frame3_spec <- compare_frame3_spec %>% filter(!AccID %in% nChars2)
    
    #no adjustment
    species_performance <- performance_single_run(compare_frame3_spec, level="Species")
    species_performance <- species_performance[-which(rownames(species_performance)=="Unassigned Correct"),,drop=FALSE] 
    #fix it
    rownames(species_performance)[which(rownames(species_performance)=="Unassigned Incorrect")] <- "Unassigned"
    
    write.table(species_performance, file= file.path(outputDir, paste0(regionName, "_Species_performance.tsv")), sep="\t", col.names = NA)
    save_file <- list("Spec"=compare_frame3_spec, "Genus"=compare_frame3)
    save(save_file, file= file.path(outputDir, paste0(regionName, "_full_comparisons.RData")))
    
  }
  
  #now calcualte the adjusted stats
  if(!historic){
    ##read in seedData
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
    ## would need to update this line since accessions are not consistent between GTDB releases 
    ## if we want to bench historically.
    taxdata_seed <- taxdata %>% filter(primaryAccession %in% SeedTax$primaryAccession)
    taxdata_SP <- word(taxdata_seed$Species, 1, 2)
    
    if(SILVA){
      #list of species that are in the seed database
      taxdata_SP <- taxdata_SP[-which(is.na(taxdata_SP))]
    }
    
    #need to do this for both species and genus
    compare_frame3 <- compare_frame3 %>% mutate(seed_genus=Genus %in% taxdata_seed$Genus)
    compare_frame3_spec <- compare_frame3_spec %>% mutate(seed_species=Species %in% taxdata_SP)
    
    #Make a corrected genus flag where we set previous set to NA. Set to true if they are unassigned and not in the seed DB
    compare_frame3 <- compare_frame3 %>% mutate(Flag.genus_cor=ifelse(is.na(Flag.genus) & !seed_genus, TRUE, Flag.genus))
    
    #Same as above for species
    compare_frame3_spec <- compare_frame3_spec %>% mutate(Flag_cor=ifelse(is.na(Flag) & !seed_species, TRUE, Flag))

    #switch so that we can use in performance table function
    compare_frame3$Flag.genus <- compare_frame3$Flag.genus_cor
    
    compare_frame3_spec$Flag <- compare_frame3_spec$Flag_cor
    
    #get adjusted species performance
    t3 <- performance_single_run(compare_frame3_spec, "Species")
    write.table(t3, file= file.path(outputDir, paste0(regionName, "_Species_performance_adjust.tsv")), sep="\t", col.names = NA)
    #get adjusted genus performance
    t4 <- performance_single_run(compare_frame3, "Genus")
    write.table(t4, file= file.path(outputDir, paste0(regionName, "_Genus_performance_adjust.tsv")), sep="\t", col.names = NA)
    
    save_file <- list("Spec"=compare_frame3_spec, "Genus"=compare_frame3)
    save(save_file, file= file.path(outputDir, paste0(regionName, "_full_comparisons_adjust.RData")))
    
     
    
  }
  
  
}


performance_single_run <- function(compareData, level){
  
  if(level=="Genus"){
    TP <- compareData %>% filter(Flag.genus) %>% nrow()
    FP <- compareData %>% filter(!Flag.genus) %>% nrow()
    TN <- 0
    FN <- compareData %>% filter(is.na(Genus_Assignment) & is.na(Flag.genus)) %>% nrow()
    
    multCorrect <- compareData %>% filter(!Genus_Assignment==Genus & Flag.genus) %>% nrow() /nrow(compareData)
    uniqueCorrect <- compareData %>% filter(Genus_Assignment==Genus) %>% nrow() /nrow(compareData)
    unassignedCorrect <- compareData %>% filter(is.na(Genus_Assignment) & Flag.genus) %>% nrow() /nrow(compareData)
    unassignedIncorrect <- compareData %>% filter(is.na(Genus_Assignment) & is.na(Flag.genus)) %>% nrow() /nrow(compareData)
  }
  
  if(level=="Species"){
    TP <- compareData %>% filter(Flag) %>% nrow()
    FP <- compareData %>% filter(!Flag) %>% nrow()
    TN <- 0
    FN <- compareData %>% filter(is.na(Species_Assignment) & is.na(Flag)) %>% nrow()
    
    multCorrect <- compareData %>% filter(!Species_Assignment==Species & Flag) %>% nrow() /nrow(compareData)
    uniqueCorrect <- compareData %>% filter(Species_Assignment==Species) %>% nrow() /nrow(compareData)
    unassignedCorrect <- compareData %>% filter(is.na(Species_Assignment) & Flag) %>% nrow() /nrow(compareData)
    unassignedIncorrect <- compareData %>% filter(is.na(Species_Assignment) & is.na(Flag)) %>% nrow() /nrow(compareData)
  }
  
  #Calculate metrics
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  f1 <- 2 * (precision * recall) / (precision + recall)
  fpr <- FP / nrow(compareData)
  
  
  #set up the output table
  rows1 <- c(             "Accuracy", "Precision", "Recall", "F1 Score",
                          "Uniquely Correct", "One-to-many Correct", "Incorrect", "Unassigned Correct", "Unassigned Incorrect")
  
  table.out <- matrix(NA, nrow = length(rows1), ncol=1)
  colnames(table.out) <- c("Query_Tool")
  rownames(table.out) <- rows1
  table.out["Uniquely Correct", "Query_Tool"] <- uniqueCorrect
  table.out["One-to-many Correct", "Query_Tool"] <- multCorrect
  table.out["Incorrect", "Query_Tool"] <- fpr
  table.out["Unassigned Correct", "Query_Tool"] <- unassignedCorrect
  table.out["Unassigned Incorrect", "Query_Tool"] <- unassignedIncorrect
  
  table.out["Accuracy", "Query_Tool"] <- accuracy
  table.out["Precision", "Query_Tool"] <- precision
  table.out["Recall", "Query_Tool"] <- recall
  table.out["F1 Score", "Query_Tool"] <- f1
  
  print(round(table.out, 3))
}
