### Generate query sequences that are in genera that are not in the seed database
require(docopt)

'Usage:
  Generate_novel_genera_query_dev.R [-p <parathaa_PATH> -s <seed_taxonomy> -t <silva_taxonomy_file>]
  
  
Options:
  -p directory where parathaa github repo is cloned
  -s seed_db_tax
  -t assignments
  ]' -> doc




opts <- docopt(doc)



library(tidyr)
library(dplyr)
library(seqinr)

parathaaDir <- (opts$p)
source(file.path(parathaaDir, "parathaa/utility/SILVA.species.editor.dev.R"))

## read in the seed database
taxonomyRanks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
seedTaxFile <- opts$s

SeedTax <- read.table(seedTaxFile , header=F, fill=TRUE,sep='\t')

suppressWarnings({
  SeedTax <- SeedTax %>%
    separate(col=V1, into=c("primaryAccession", "ArbID"), sep="\\.") %>%
    separate(col=V2, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep=";") %>%
    filter(Kingdom=="Bacteria" & !is.na(Genus) & Genus!="")
})



#Seed database contains 2076 genera

#load in full silva database

fullSilvaTaxFile <- opts$t


#Load in the taxdata in the format that we want to use for benchmarking to remove unnamed and metagenome queries
taxdata <- read.table(fullSilvaTaxFile , header=T, fill=TRUE,sep='\t', quote="")
taxdata <- taxdata %>%
  unite("AccID", c("primaryAccession", "start", "stop"), sep=".", remove=F)
taxdata <- taxdata %>%
  mutate(taxonomy=paste0(path, organism_name))

suppressWarnings({
  taxdata <- taxdata %>%
    select(AccID, primaryAccession, start, stop, taxonomy) %>%
    separate(col=taxonomy, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
})


## Additional changes (getting rid of subspecies)
taxdata <- SILVA.species.editor(taxdata)

#filter it to the taxa that are found in see
taxdata_seed <- taxdata %>% filter(primaryAccession %in% SeedTax$primaryAccession)

#this gets all the genera we cannot include
Seed_Genera <- unique(taxdata_seed$Genus)



non_seed_Genera <- taxdata %>% filter(!Genus %in% Seed_Genera)


### random select 10,000 sequences
set.seed(1995)
new_query_ids <- sample(non_seed_Genera$AccID, size = 10000)


write.table(new_query_ids, file = "input/novel_query_ids.txt", col.names = F, row.names = F, quote=F)


