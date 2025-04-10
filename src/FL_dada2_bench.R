#Script that runs FL benchmark
require(docopt)

'Usage:
  FL_dada2_bench.R [-p <parathaa_PATH> --dada_db_FL <dada2_db> -t <input_taxonomy> -o <output> --paraAssign <para_Assignments> --query <query seqs> -s <seed_data> -b <min boot>]
  
  
Options:
  -p directory where parathaa github repo is cloned
  --dada_db_FL location of DADA2 db
  --paraAssign location of parathaa assignments
  --query location of query reads 
  -t input taxonomy
  -o output
  -s seed_data
  -b minboot [default=80]
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


source("src/performance.table.R")

message(paste0("running dada2 full length benchmark with minboot: ", as.numeric(opts$b)))

run.synthetic.data(parathaaFile = opts$paraAssign, 
               sequenceFile = opts$query, 
               outputDir = opts$o, 
               DADAdb = opts$dada_db_FL, 
               inFileTaxdata = opts$t, 
               inFileSeedDB = opts$s,
               minboot = as.numeric(opts$b), regionName = "FL", SILVA = T, full_length = T, historic = F)