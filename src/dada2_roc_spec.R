### ROC Curve generation

## We will use the dada2 FL database as we want to only use naive bayes on all of these as exact
## species matching does not allow for ROC curve creation as there is no parameter to adjust


### First need to generate dada2 classifications with 0 minBoot...
require(docopt)

'Usage:
    dada2_roc_spec.R [-p <parathaa_PATH>  --dada_db_FL <dada2_db_FL> -t <input_taxonomy> -o <output> --queryV1V2 <query_seqs> --queryV4V5 <query_seqs> --queryFL <query_seqs> ]
  
  
Options:
  --queryV1V2 V1V2 querys
  --queryV4V5 V4V5 querys
  --queryFL full length querys
  --dada_db_FL location of dada FL db
  -o output output directory
  ]' -> doc


opts <- docopt(doc)

parathaaDir <- (opts$p)

library(dada2)
library(seqinr)

classify_dada2_seqs <- function(seqs, regionName){
  set.seed(1995)
  taxa <- assignTaxonomy(seqs, opts$dada_db_FL, minBoot = 0, outputBootstraps = T, multithread = T)
  
  ## First, get names and sequences from fasta file
  getNames <- read.fasta(file(seqs), as.string = TRUE,
                         forceDNAtolower = FALSE, whole.header = FALSE)
  #get the names and split it by (tab as we only want to keep the accession ID)
  names1 <- str_split(getName(getNames), "\t", simplify=TRUE)
  #remove > in fasta headers if they exist
  names1 <- names1[,1] %>%
    str_remove(">")
  #create data frame with sequence and ID
  name.df <- data.frame("sequence" = unlist(getSequence(getNames, as.string=T)), taxaIDs = names1)
  
  stopifnot(length(which(name.df$sequence == rownames(taxa$tax))) == dim(taxa$tax)[1])
  taxa$names <- name.df$taxaIDs
  saveRDS(taxa, paste0(opts$o, regionName, "_dada2_class_boot0.RDS"))
}

classify_dada2_seqs(opts$queryV1V2, regionName = "V1V2")
classify_dada2_seqs(opts$queryV4V5, regionName="V4V5")
classify_dada2_seqs(opts$queryFL, regionName="FL")
