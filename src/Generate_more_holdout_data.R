
require(docopt)

'Usage:
  Generate_more_holdout_data.R [-p <parathaa_PATH> -s <seed_taxonomy> -t <silva_taxonomy_file> -q <old_queries>]
  
  
Options:
  -p directory where parathaa github repo is cloned
  -s seed_db_tax
  -t assignments
  -q old_queries
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

suppressWarnings({
  seedTax <- seedTax %>% separate(col=V2, into=taxonomyRanks, sep = ";")
})

seedTax <- seedTax %>% mutate(primaryaccession=gsub("\\..*", "", V1))

#Seed database contains 2076 genera

#load in full silva database

fullSilvaTaxFile <- opts$t

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
taxdata <- taxdata %>% filter(!is.na(Species))


#filter out the sequences that in the seed database
taxdata_filt <- taxdata %>% filter(!primaryAccession %in% seedTax$primaryaccession)


### filter out the ids in the original benchmark
prev_test_reads <- read.table(file=opts$q)

#filter them out
taxdata_filt_filt <- taxdata_filt %>% filter(!AccID %in% prev_test_reads$V1)


#remove euks
taxdata_filt_filt <- taxdata_filt_filt %>% filter(Kingdom!="Eukaryota")

### sample 30,000 sequences then generate three datasets
set.seed(1995)
hold_out_querys <- sample(taxdata_filt_filt$AccID, size=30000, replace=F)

hold_outs_1 <- hold_out_querys[1:10000]
hold_outs_2 <- hold_out_querys[10001:20000]
hold_outs_3 <- hold_out_querys[20001:30000]


write.table(hold_outs_1, file = "input/holdout_query1.txt", col.names = F, row.names = F, quote=F)
write.table(hold_outs_2, file = "input/holdout_query2.txt", col.names = F, row.names = F, quote=F)
write.table(hold_outs_3, file = "input/holdout_query3.txt", col.names = F, row.names = F, quote=F)

### non-random holdout testing replicate to 

#taxdata seedless are taxa that are not in the seed database are not in the original benchmark
#but have genus that are contained within the seed database (just like how the original benchmark was done)
#only issue with this is that in some cases we will have less sequences from genera that have a low number of sequences within them 
#(i.e. missing those that have < 5 reps)

taxdata_seedless <- taxdata_filt_filt %>% 
  filter(!primaryAccession %in% seedTax$primaryAccession) %>%
  ## Subset to Genera in Seed for now
  filter(Genus %in% unique(seedTax$Genus))

GenusCounts <- taxdata_seedless %>% dplyr::count(Genus)
GenusCounts$subsetN <- pmin(pmax(20, ceiling(GenusCounts$n*0.01)), GenusCounts$n)
## Include only genera from seed for now
GenusCounts <- GenusCounts %>% filter(Genus %in% unique(seedTax$Genus))


library(purrr)
set.seed(978)
subsample <- taxdata_seedless %>% 
  group_split(Genus) %>% 
  map2_dfr(GenusCounts$subsetN, ~ slice_sample(.x, n = .y))

## Print IDs to file
readr::write_lines(subsample$AccID, file="input/original_holdout1.txt")

