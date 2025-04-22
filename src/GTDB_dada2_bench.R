#### R script to run the GTDB benchmarks

require(docopt)

'Usage:
    GTDB_dada2_bench.R [-p <parathaa_PATH> --dada_db <dada2_db> --dada_db_sp <dada2_db_sp> --dada_db_FL <dada2_db_FL> -t <input_taxonomy> -o <output> --paraAssignV1V2 <taxonomy_file> --paraAssignV4V5 <taxonomy_file> --paraAssignFL <taxonomy_file> --queryV1V2 <query seqs> --queryV4V5 <query seqs> --queryFL <query seqs>]
  
  
Options:
  -p directory where parathaa github repo is cloned
  --paraAssignV1V2 
  --paraAssignV4V5
  --paraAssignFL 
  --queryV1V2
  --queryV4V5
  --queryFL
  --dada_db location of species classifer
  --dada_db_sp number of threads
  --dada_db_FL
  -t input taxonomy
  -o output
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

##There is probably a better way of sourcing this file... but for now leave it as is...
source("src/performance.table.R")


DADAdb <- opts$dada_db
DADAdb.sp <- opts$dada_db_sp
inFileTaxdata <- opts$t


##### CALL FUNCTION ######

## V1V2
message('Running V1V2 bench')
run.synthetic.data(parathaaFile = as.character(opts$paraAssignV1V2), 
                   sequenceFile = opts$queryV1V2,
                   regionName = "V1V2", 
                   outputDir=paste0(opts$o, "/V1V2/", sep=""),
                   dadaAllowMult = T,
                   DADAdb = DADAdb,
                   DADAdb.sp = DADAdb.sp,
                   inFileTaxdata = inFileTaxdata, SILVA = F, historic = T)


# V4V5
message("Running V4V5 bench")
run.synthetic.data(parathaaFile = opts$paraAssignV4V5, 
                   sequenceFile = opts$queryV4V5,
                   regionName = "V4V5", 
                   outputDir=paste0(opts$o, "/V4V5/", sep=""),
                   dadaAllowMult = T,
                   DADAdb = DADAdb,
                   DADAdb.sp = DADAdb.sp,
                   inFileTaxdata = inFileTaxdata, SILVA = F, historic = T)

# Full length
message("Running FL bench")
run.synthetic.data(parathaaFile = opts$paraAssignFL, 
                   sequenceFile = opts$queryFL,
                   regionName = "FL", 
                   outputDir=paste0(opts$o, "/FL/", sep=""),
                   dadaAllowMult = T,
                   DADAdb = opts$dada_db_FL,
                   inFileTaxdata = inFileTaxdata, SILVA = F, historic = T, full_length = T, minboot = 80)

# Full length exact
message("Running FL exact Bench")
run.synthetic.data(parathaaFile = opts$paraAssignFL, 
                   sequenceFile = opts$queryFL,
                   regionName = "FL", 
                   outputDir=paste0(opts$o, "/FL_exact/", sep=""),
                   dadaAllowMult = T,
                   DADAdb = DADAdb,
                   DADAdb.sp = DADAdb.sp,
                   inFileTaxdata = inFileTaxdata, SILVA = F, historic = T, full_length = F)
