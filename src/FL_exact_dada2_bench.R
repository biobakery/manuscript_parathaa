### Full length exact matching
require(docopt)

'Usage:
    FL_exact_dada2_bench.R [-p <parathaa_PATH> --dada_db <dada2_db> --dada_db_sp <dada2_db_sp> -t <input_taxonomy> -o <output> --paraAssign <taxonomy_file> --query <query seqs> -s <seed_data>]
  
  
Options:
  -p directory where parathaa github repo is cloned
  --dada_db location of DADA2 db
  --dada_db_sp location of DADA2 species db
  --paraAssign parathaa assignments for V4V5
  --query query sequence file for V4V5
  -t input taxonomy
  -o output
  -s seed_data
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
source(file.path(parathaaDir, "parathaa/utility/SILVA.species.editor.dev.R"))

##There is probably a better way of sourcing this file... but for now leave it as is...
source("src/performance.table.R")


DADAdb <- opts$dada_db
DADAdb.sp <- opts$dada_db_sp
inFileTaxdata <- opts$t
inFileSeedDB <- opts$s


##### CALL FUNCTION ######
run.synthetic.data(parathaaFile = opts$paraAssign, 
                   sequenceFile = opts$query,
                   regionName = "FL", 
                   outputDir=paste(opts$o, "/Figures/synth_mult_arc/", sep=""),
                   dadaAllowMult = T,
                   DADAdb = DADAdb,
                   DADAdb.sp = DADAdb.sp,
                   inFileTaxdata = inFileTaxdata,
                   inFileSeedDB = inFileSeedDB, SILVA = T)

run.synthetic.data(parathaaFile = opts$paraAssign, 
                   sequenceFile = opts$query,
                   regionName = "FL", 
                   outputDir=paste(opts$o, "/Figures/synth_nomult_arc/", sep=""),
                   dadaAllowMult = F,
                   DADAdb = DADAdb,
                   DADAdb.sp = DADAdb.sp,
                   inFileTaxdata = inFileTaxdata,
                   inFileSeedDB = inFileSeedDB, SILVA = T)
