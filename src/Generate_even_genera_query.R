### Generate query sequences that are in genera that are not in the seed database

require(docopt)

'Usage:
  Generate_even_genera_query.R [-p <parathaa_PATH> -s <seed_taxonomy> -t <silva_taxonomy_file>]
  
  
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


seedTax <- read.table(seedTaxFile, header=F, sep="\t")


seedTax <- seedTax %>% separate(col=V2, into=taxonomyRanks, sep = ";")
seedTax <- seedTax %>% mutate(primaryaccession=gsub("\\..*", "", V1))

SeedGenera <- unique(seedTax$Genus)

#Seed database contains 2076 genera

#load in full silva database

fullSilvaTaxFile <- opts$t

taxdata <- read.table(fullSilvaTaxFile , header=T, fill=TRUE,sep='\t', quote="")
taxdata <- taxdata %>%
  unite("AccID", c("primaryAccession", "start", "stop"), sep=".", remove=F)
taxdata <- taxdata %>%
  mutate(taxonomy=paste0(path, organism_name))

taxdata <- taxdata %>%
  select(AccID, primaryAccession, start, stop, taxonomy) %>%
  separate(col=taxonomy, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")

## Additional changes (getting rid of subspecies)
taxdata <- SILVA.species.editor(taxdata)
taxdata <- taxdata %>% filter(!is.na(Species))

### Sample so there is only 1 query per genera and that sequence is not in the seed database
taxdata_filt <- taxdata %>% filter(!primaryAccession %in% seedTax$primaryaccession)

set.seed(1995)
genera_sampling <- taxdata_filt %>% group_by(Genus) %>% sample_n(size=1) %>% ungroup()

new_query_ids <- genera_sampling$AccID


write.table(new_query_ids, file = "input/even_queryIDs.txt", col.names = F, row.names = F, quote=F)


