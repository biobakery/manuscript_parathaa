
require(docopt)

'Usage:
  Generate_more_holdout_data.R [-p <parathaa_PATH> -s <seed_taxonomy> -t <silva_taxonomy_file> -q <old_queries>]
  
  
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

taxonomyRanks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
seedTaxFile <- opts$s


seedTax <- read.table(seedTaxFile, header=F, sep="\t")


seedTax <- seedTax %>% separate(col=V2, into=taxonomyRanks, sep = ";")
seedTax <- seedTax %>% mutate(primaryaccession=gsub("\\..*", "", V1))

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


#filter out the sequences that in the seed database
taxdata_filt <- taxdata %>% filter(!primaryAccession %in% seedTax$primaryaccession)


### filter out the ids in the original benchmark
prev_test_reads <- read.table(opts$q)

#filter them out
taxdata_filt_filt <- taxdata_filt %>% filter(!AccID %in% prev_test_reads$V1)


### sample 20,000 sequences then generate two datasets
set.seed(1995)
hold_out_querys <- sample(taxdata_filt_filt$AccID, size=20000, replace=F)

hold_outs_1 <- hold_out_querys[1:10000]
hold_outs_2 <- hold_out_querys[10001:20000]


write.table(hold_outs_1, file = "input/holdout_query1.txt", col.names = F, row.names = F, quote=F)
write.table(hold_outs_2, file = "input/holdout_query2.txt", col.names = F, row.names = F, quote=F)

