####################################
# Create Synthetic reads for benchmarking
####################################


require(docopt)

'Usage:
  Synthetic.Reads.R [-p <parathaa_PATH> -s <seed_taxonomy> -t <silva_taxonomy_file>]
  
  
Options:
  -p directory where parathaa github repo is cloned
  -s seed_db_tax
  -t assignments
  ]' -> doc




opts <- docopt(doc)


parathaaDir <- (opts$p)

source(file.path(parathaaDir, "parathaa/utility/SILVA.species.editor.dev.R"))

## Read in SILVA 138.1 taxonomy for subsetting
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
inFileSeedDB <- opts$s
SeedTax <- read.table(inFileSeedDB , header=F, fill=TRUE,sep='\t')

SeedTax <- SeedTax %>%
  separate(col=V1, into=c("primaryAccession", "ArbID"), sep="\\.") %>%
  separate(col=V2, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep=";") %>%
  filter(Kingdom=="Bacteria" & !is.na(Genus) & Genus!="") 

taxdata_seedless <- taxdata %>% 
  filter(!primaryAccession %in% SeedTax$primaryAccession) %>%
  ## Subset to Genera in Seed for now
  filter(Genus %in% unique(SeedTax$Genus))


## How many to sample from each Genus
GenusCounts <- taxdata_seedless %>% dplyr::count(Genus)
GenusCounts$subsetN <- pmin(pmax(20, ceiling(GenusCounts$n*0.01)), GenusCounts$n)
## Include only genera from seed for now
GenusCounts <- GenusCounts %>% filter(Genus %in% unique(SeedTax$Genus))


library(purrr)
set.seed(978)
subsample <- taxdata_seedless %>% 
  group_split(Genus) %>% 
  map2_dfr(GenusCounts$subsetN, ~ slice_sample(.x, n = .y))

## Print IDs to file
readr::write_lines(subsample$AccID, file="input/subsampleIDs_SeedGenera.txt")
