#### Script that runs V1V2, V4V5 benchmarks
require(docopt)

'Usage:
  V1V2_V4V5_dada2_bench.R [-p <parathaa_PATH> --dada_db <dada2_db> --dada_db_sp <dada2_db_sp> -t <input_taxonomy> -o <output> --paraAssignV4V5 <V4V5_taxonomy_file> --paraAssignV1V2 <V1V2_taxonomy_file> --queryV4V5 <V4V5 seqs> --queryV1V2 <V1V2 seqs> -s <seed_data>]
  
  
Options:
  -p directory where parathaa github repo is cloned
  --dada_db location of DADA2 db
  --dada_db_sp location of DADA2 species db
  --paraAssignV4V5 parathaa assignments for V4V5
  --paraAssignV1V2 parathaa assignments for V1V2
  --queryV4V5 query sequence file for V4V5
  --queryV1V2 query sequence file for V1V2
  -t input taxonomy
  -o output
  -s seed_data
  ]' -> doc


opts <- docopt(doc)

parathaaDir <- (opts$p)

library(phyloseq) 
library(dplyr) 
library(stringr)
suppressPackageStartupMessages(library(ggtree))
library(treeio)
library(dada2)
library(ggplot2)
library(castor)
suppressPackageStartupMessages(library(phytools))
#library(flextable)
suppressPackageStartupMessages(library(seqinr))
source(file.path(parathaaDir, "parathaa/utility/SILVA.species.editor.dev.R"))

##There is probably a better way of sourcing this file... but for now leave it as is...
source("src/performance.table.R")


## Define variables used across all function calls
## 
DADAdb <- opts$dada_db
DADAdb.sp <- opts$dada_db_sp
inFileTaxdata <- opts$t
inFileSeedDB <- opts$s


##### CALL FUNCTION ######
run.synthetic.data(parathaaFile = opts$paraAssignV4V5, 
                   sequenceFile = opts$queryV4V5,
                   regionName = "V4V5", 
                   outputDir=paste(opts$o, "/Figures/synth_mult_arc/", sep=""),
                   dadaAllowMult = T,
                   DADAdb = DADAdb,
                   DADAdb.sp = DADAdb.sp,
                   inFileTaxdata = inFileTaxdata,
                   inFileSeedDB = inFileSeedDB, SILVA = T, historic = F, full_length = F, minboot = 50)

run.synthetic.data(parathaaFile = opts$paraAssignV1V2, 
                   sequenceFile = opts$queryV1V2, 
                   regionName = "V1V2", 
                   outputDir=paste(opts$o, "/Figures/synth_mult_arc", sep=""),
                   dadaAllowMult = T,
                   DADAdb = DADAdb,
                   DADAdb.sp = DADAdb.sp,
                   inFileTaxdata = inFileTaxdata,
                   inFileSeedDB = inFileSeedDB, SILVA=T, historic=F, full_length = F, minboot=50)

run.synthetic.data(parathaaFile = opts$paraAssignV4V5, 
                   sequenceFile = opts$queryV4V5,
                   regionName = "V4V5", 
                   outputDir=paste(opts$o, "/Figures/synth_nomult_arc/", sep=""),
                   dadaAllowMult = F,
                   DADAdb = DADAdb,
                   DADAdb.sp = DADAdb.sp,
                   inFileTaxdata = inFileTaxdata,
                   inFileSeedDB = inFileSeedDB, SILVA=T, historic = F, full_length = F, minboot = 50)

run.synthetic.data(parathaaFile = opts$paraAssignV1V2, 
                   sequenceFile = opts$queryV1V2, 
                   regionName = "V1V2", 
                   outputDir=paste(opts$o, "/Figures/synth_nomult_arc", sep=""),
                   dadaAllowMult = F,
                   DADAdb = DADAdb,
                   DADAdb.sp = DADAdb.sp,
                   inFileTaxdata = inFileTaxdata,
                   inFileSeedDB = inFileSeedDB, SILVA=T, historic = F, full_length = F, minboot=50)