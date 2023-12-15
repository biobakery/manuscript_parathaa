
####################################
# Create dada2 compatible database
# Start with degapped seed db
####################################
parathaaDir <- ("C:/Users/mshort/Documents/proj/parathaa/")
source(file.path(parathaaDir, "parathaa/utility/SILVA.species.editor.dev.R"))

## Read in seed db 
seed.db <- seqinr::read.fasta("input/silva.seed_v138_1.ng.fasta", forceDNAtolower = FALSE,
                              as.string=TRUE)


#Extract names: "Accession_number.Arb_ID"
nam <- getName(seed.db)
# Extract Accession numbers
nam2 <- strsplit(nam, split="[.]")
nam3 <- sapply(nam2,"[[",1)

if(FALSE){
## Get dada2 species assignment file for larger SILVA db
silva.sp <- seqinr::read.fasta("input/silva_species_assignment_v138.1.fa" ,forceDNAtolower = FALSE,
                               as.string=TRUE)
## Extract names: "Accession.number.start.end"
silva.nam <- getName(silva.sp)
silva.nam2 <- strsplit(silva.nam, split="[.]")
# Extract Accession numbers
silva.nam3 <- sapply(silva.nam2,"[[",1)

# How many seed database sequences are found in the dada2 silva species naming file?
length(which(nam3 %in% silva.nam3))

# Which seqs are in seed db but not in silva_species_assignment training file?
excluded <- seed.db[which(!nam3 %in% silva.nam3)]
excluded <- nam3[!nam3 %in% silva.nam3]
taxdata %>% filter(primaryAccession %in% excluded) %>% View

## Short answer: Eukaryotes, and species without names (some of which have numbers, i.e. Paenibacillus sp. Cp_S316	45277)
## We want some of these (the named ones that are not Eukaryotes), so we try a different approach

### Before doing this, going to try using only those from sp naming file:
#Seqs that are in dada2's silva-based species naming file that also show up in the seed database.
silva.sp.seed <- silva.sp[silva.nam3 %in% nam3]
silva.sp.seed.names <- substring(getAnnot(silva.sp.seed), 2)

write.fasta(silva.sp.seed, names=silva.sp.seed.names, nbchar=80, file.out = "/Users/mis696/proj/parathaa/input/silva_species_assignment_v138.1.seedonly.fa")

### And similarly, take full training set and subset to seed db:
silva.train <- seqinr::read.fasta("/Users/mis696/proj/parathaa/input/silva_nr99_v138.1_train_set.fa" ,forceDNAtolower = FALSE,
                                  as.string=TRUE)

silva.train.seed <- silva.train[silva.train %in% silva.sp.seed]
write.fasta(silva.train.seed, names=getName(silva.train.seed), nbchar=80, file.out = "/Users/mis696/proj/parathaa/input/silva_nr99_v138.1_train_set.seedonly.fa")
}

# Read in taxonomy file from full SILVA db
##tax <- read.delim("input/silva.seed_v138_1.tax", header=FALSE, row.names = 1, col.names = c("Name", "Taxonomy"))

inFileTaxdata <- "input/taxmap_slv_ssu_ref_138.1.txt"


taxdata <- read.table(inFileTaxdata , header=T, fill=TRUE,sep='\t', quote="")
taxdata <- taxdata %>%
  unite("AccID", c("primaryAccession", "start", "stop"), sep=".", remove=F)
taxdata <- taxdata %>%
  separate(col=path, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep=";") %>%
  dplyr::rename(Species = organism_name) %>%
  filter(Kingdom!="Eukaryota") ## NOTE to Jacob: This was originally limiting to Bacteria


taxdata.sp <- taxdata %>% group_by(primaryAccession, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% tally()

#Subset to primaryAccessions from seed db
taxmatch <- taxdata.sp[which(taxdata.sp$primaryAccession %in% nam3),] 
taxmatch.df <- as.data.frame(taxmatch)
rownames(taxmatch.df) <- taxmatch.df$primaryAccession

taxmatch2 <- taxmatch.df[nam3,]

#Out of curiosity, find doubles in taxdata
#doubles <- names(table(taxdata.sp$primaryAccession)[which(table(taxdata.sp$primaryAccession)>1)])
#taxdata.sp[which(taxdata.sp$primaryAccession %in% doubles),] %>% View()

## More seed sequences show up in this full database than in the pre-trained database (based on NR99)
dim(taxmatch2)
dim(taxmatch2[!is.na(taxmatch2$Species),])
taxmatch3 <- SILVA.species.editor(taxmatch2)
taxmatch4 <- taxmatch3 %>% filter(!is.na(Species))

## Write seed db in training file format for DADA2 
#seed.db.species <- seed.db[which(!is.na(taxmatch2$Species))]
seed.db.species <- seed.db[which(taxmatch2$primaryAccession %in% taxmatch4$primaryAccession)]
seed.db.species.names <- paste(nam3[which(!is.na(taxmatch3$Species))], taxmatch3$Species[which(!is.na(taxmatch3$Species))])

#write.fasta(sequences = seed.db.species, names=seed.db.species.names, nbchar=80, file.out = "input/20231130_silva.seed_v138_1.ng.dada.sp.fasta")
write.fasta(sequences = seed.db.species, names=seed.db.species.names, nbchar=80, file.out = "input/20231215_silva.seed_v138_1.ng.dada.sp.fasta")


### Get full taxonomies for seed database (Bacteria only)
bacRows <- which(!is.na(taxmatch2$Kingdom))
seed.db.bac <- seed.db[bacRows]
seed.db.bac.names <- apply(cbind(taxmatch2$Kingdom[bacRows],
                                 taxmatch2$Phylum[bacRows],
                                 taxmatch2$Class[bacRows],
                                 taxmatch2$Order[bacRows],
                                 taxmatch2$Family[bacRows],
                                 taxmatch2$Genus[bacRows]), 1, 
                           function(x) paste(x[!is.na(x)], collapse = ";"))
seed.db.bac.names <- paste0(seed.db.bac.names, ";")


#write.fasta(sequences = seed.db.bac , names = seed.db.bac.names, nbchar = 80, file.out = "input/20230111.silva.seed_v138_1.ng.dada.fasta")
write.fasta(sequences = seed.db.bac , names = seed.db.bac.names, nbchar = 80, file.out = "input/20231215.silva.seed_v138_1.ng.dada.fasta")
