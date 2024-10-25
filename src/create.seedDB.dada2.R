
####################################
# Create dada2 compatible database
# Start with degapped seed db
####################################


require(docopt)

'Usage:
  create.seedDB.dada2.R [-p <parathaa_PATH> -s <seed_db> -t <silva_taxonomy_file>]
  
  
Options:
  -p directory where parathaa github repo is cloned
  -s seed_db
  -t assignments
  ]' -> doc




opts <- docopt(doc)




parathaaDir <- opts$p
source(file.path(parathaaDir, "parathaa/utility/SILVA.species.editor.dev.R"))

## Read in seed db 
seed.db <- seqinr::read.fasta(opts$s, forceDNAtolower = FALSE,
                              as.string=TRUE)


#Extract names: "Accession_number.Arb_ID"
nam <- seqinr::getName(seed.db)
# Extract Accession numbers
nam2 <- strsplit(nam, split="[.]")
nam3 <- sapply(nam2,"[[",1)

# Read in taxonomy file from full SILVA db
inFileTaxdata <- opts$t


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



## fix taxonomic labels
taxmatch3 <- SILVA.species.editor(taxmatch2)
#only grab those with species names for the species level dada alignment DB
taxmatch4 <- taxmatch3 %>% filter(!is.na(Species))

## Write seed db in training file format for DADA2 
seed.db.species <- seed.db[which(taxmatch2$primaryAccession %in% taxmatch4$primaryAccession)]
seed.db.species.names <- paste(nam3[which(!is.na(taxmatch3$Species))], taxmatch3$Species[which(!is.na(taxmatch3$Species))])

seqinr::write.fasta(sequences = seed.db.species, names=seed.db.species.names, nbchar=80, file.out = "input/20231215_silva.seed_v138_1.ng.dada.sp.fasta")


### Get full taxonomies for seed database (non-eurkaytoa)
bacRows <- which(!is.na(taxmatch2$Kingdom))
seed.db.bac <- seed.db[bacRows]
seed.db.bac.names <- apply(cbind(taxmatch3$Kingdom[bacRows],
                                 taxmatch3$Phylum[bacRows],
                                 taxmatch3$Class[bacRows],
                                 taxmatch3$Order[bacRows],
                                 taxmatch3$Family[bacRows],
                                 taxmatch3$Genus[bacRows]), 1, 
                           function(x) paste(x[!is.na(x)], collapse = ";"))
seed.db.bac.names <- paste0(seed.db.bac.names, ";")


seqinr::write.fasta(sequences = seed.db.bac , names = seed.db.bac.names, nbchar = 80, file.out = "input/20231215.silva.seed_v138_1.ng.dada.fasta")


### Get full taxonomies + Species for the full length assignment database
bacRows <- which(!is.na(taxmatch2$Kingdom))
seed.db.bac <- seed.db[bacRows]
seed.db.bac.names <- apply(cbind(taxmatch3$Kingdom[bacRows],
                                 taxmatch3$Phylum[bacRows],
                                 taxmatch3$Class[bacRows],
                                 taxmatch3$Order[bacRows],
                                 taxmatch3$Family[bacRows],
                                 taxmatch3$Genus[bacRows],
                                 taxmatch3$Species[bacRows]), 1, 
                           function(x) paste(x[!is.na(x)], collapse = ";"))
seed.db.bac.names <- paste0(seed.db.bac.names, ";")

seqinr::write.fasta(sequences = seed.db.bac , names = seed.db.bac.names, nbchar = 80, file.out = "input/20231215.silva.seed_v138_1.ng.dada_FL.fasta")
