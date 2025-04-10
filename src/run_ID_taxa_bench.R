### Script for running IDTaxa

require(docopt)

'Usage:
  run_ID_taxa_bench.R [--IDTAXA_spec_db <IDTAXA species db> --IDTAXA_genus_db <IDTAXA genus db> --queryV4V5 <V4V5 query seqs> --queryV1V2 <V1V2 query seqs> --queryFL <FL query seqs> --runFL <FALSE> -t <taxonomy file> --threads <threads> -s <seed_db> -o <output_dir> -p <parathaa install>]
  
  
Options:
  -p parathaa directory
  --IDTAXA_spec_db location of DADA2 db
  --IDTAXA_genus_db location of DADA2 species db
  --queryV4V5 query sequence file for V4V5
  --queryV1V2 query sequence file for V1V2
  --queryFL query sequence file for FL if runFL is true [default: ""]
  --runFL arugement for whether to run FL or not [default: FALSE]
  -t input taxonomy
  -o output
  -s seed_data
  --threads [default: 8]
  ]' -> doc

opts <- docopt(doc)
parathaaDir <- opts$p
source("src/performance.table.R")
source(file.path(parathaaDir, "parathaa/utility/SILVA.species.editor.dev.R"))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DECIPHER))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(purrr))

#V1V2 run
Bench_IDTAXA(sequenceFile = opts$queryV1V2, genus_class = opts$IDTAXA_genus_db, 
             species_class = opts$IDTAXA_spec_db, threads = as.numeric(opts$threads), inFileTaxdata = opts$t, 
             inFileSeedDB = opts$s, SILVA = T, outputDir = opts$o, historic = FALSE, regionName="V1V2")

#V4V5 run
Bench_IDTAXA(sequenceFile = opts$queryV4V5, genus_class = opts$IDTAXA_genus_db, 
             species_class = opts$IDTAXA_spec_db, threads = as.numeric(opts$threads), inFileTaxdata = opts$t, 
             inFileSeedDB = opts$s, SILVA = T, outputDir = opts$o, historic = FALSE, regionName="V4V5")

if(opts$runFL){
  Bench_IDTAXA(sequenceFile = opts$queryFL, genus_class = opts$IDTAXA_genus_db, 
               species_class = opts$IDTAXA_spec_db, threads = as.numeric(opts$threads), inFileTaxdata = opts$t, 
               inFileSeedDB = opts$s, SILVA = T, outputDir = opts$o, historic = FALSE, regionName="FL")
}


