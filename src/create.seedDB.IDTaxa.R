### Script to create the IDTaxa species and Genus databases
### Written by Jacob T. Nearing using code from Thomas Kuntz
require(docopt)

'Usage:
  create.seedDB.IDTaxa.R [-s <seed_db> -o <output_dir>]
  
  
Options:
  -s seed_db_FL
  -o outputDir
  ]' -> doc

opts <- docopt(doc)


library(DECIPHER)
library(stringr)
library(readr)

#Read in FL file that was created for dada2
sequenceFile <- opts$s  # 
dna <- readDNAStringSet(sequenceFile)
names(dna) <- paste("Root", names(dna), sep = ";")

class_levels <- str_split(names(dna), ";")

## remove members that are empty strings
class_levels_clean <- lapply(class_levels, function(sublist) { sublist[sublist != ""]})

seqs <- dna[lapply(class_levels_clean, length)==8]

num_removed <- length(dna) - length(seqs)

message(paste0("Removed ", num_removed, " sequences that did not have a species level classification ", 
               length(seqs), " sequences remain"))

seqs <- RemoveGaps(seqs)  # there aren't any gaps, but whatever, good to make sure
seqs <- OrientNucleotides(seqs)  # also should be fine to skip

# obtain the taxonomic assignments
groups <- names(seqs) # sequence names
# assume the taxonomy begins with 'Root;'
groups <- gsub("(.*)(Root;)", "\\2", groups) # extract the group label
groupCounts <- table(groups)
u_groups <- names(groupCounts) # unique groups

message(paste0("There are ", length(u_groups), " taxonomic groups"))

# prune
maxGroupSize <- 10 # max sequences per label (>= 1)
remove <- logical(length(seqs))
for (i in which(groupCounts > maxGroupSize)) {
  index <- which(groups==u_groups[i])
  keep <- sample(length(index),
                 maxGroupSize)
  remove[index[-keep]] <- TRUE
}
message(paste0("Pruned ", sum(remove), " sequences that were over represented as suggested in IDTaxa tutorial"))

maxIterations <- 10 # must be >= 1; didn't quite converge at 10 interations, but 3 problem sequences is very few so idk maybe try more?
allowGroupRemoval <- FALSE
probSeqsPrev <- integer() # suspected problem sequences from prior iteration
for (i in seq_len(maxIterations)) {
  cat("Training iteration: ", i, "\n", sep="")
  # train the classifier
  trainingSet <- LearnTaxa(seqs[!remove],
                           names(seqs)[!remove])  # not bothering to make the taxid file, i'm still unsure what it even does (and it is optional)
  # look for problem sequences
  probSeqs <- trainingSet$problemSequences$Index
  if (length(probSeqs)==0) {
    cat("No problem sequences remaining.\n")
    break
  } else if (length(probSeqs)==length(probSeqsPrev) &&
             all(probSeqsPrev==probSeqs)) {
    cat("Iterations converged.\n")
    break
  }
  if (i==maxIterations)
    break
  probSeqsPrev <- probSeqs
  # remove any problem sequences
  index <- which(!remove)[probSeqs]
  remove[index] <- TRUE # remove all problem sequences
  if (!allowGroupRemoval) {
    # replace any removed groups
    missing <- !(u_groups %in% groups[!remove])
    missing <- u_groups[missing]
    if (length(missing) > 0) {
      index <- index[groups[index] %in% missing]
      remove[index] <- FALSE # don't remove
    }
  }
  sum(remove) # total number of sequences eliminated
}

message("A total of ", length(probSeqs), " remain after 10 rounds of iterative training")


message("Saving Species level IDTaxa database")
trainingSet %>% write_rds(paste0(opts$o,"/IdTaxa_20231215.silva.seed_v138_1.ng_FL_sp.RData"))  # again, make sure this is right



### Now write the genus level one
seqs <- dna[lapply(class_levels_clean, length)>=7]

#remove species level and only keep genus 
names(seqs) <- word(names(seqs), 1, 7, sep = ";")

num_removed <- length(dna) - length(seqs)

message(paste0("Removed ", num_removed, " sequences that did not have a genus level classification ", 
               length(seqs), " sequences remain"))

seqs <- RemoveGaps(seqs)  # there aren't any gaps, but whatever, good to make sure
seqs <- OrientNucleotides(seqs)  # also should be fine to skip

# obtain the taxonomic assignments
groups <- names(seqs) # sequence names
# assume the taxonomy begins with 'Root;'
groups <- gsub("(.*)(Root;)", "\\2", groups) # extract the group label
groupCounts <- table(groups)
u_groups <- names(groupCounts) # unique groups

message(paste0("There are ", length(u_groups), " unique taxonomic groups"))

# prune
maxGroupSize <- 10 # max sequences per label (>= 1)
remove <- logical(length(seqs))
for (i in which(groupCounts > maxGroupSize)) {
  index <- which(groups==u_groups[i])
  keep <- sample(length(index),
                 maxGroupSize)
  remove[index[-keep]] <- TRUE
}
message(paste0("Pruned ", sum(remove), " sequences that were over represented as suggested in IDTaxa tutorial"))


maxIterations <- 10 # must be >= 1; didn't quite converge at 10 interations, but 3 problem sequences is very few so idk maybe try more?
allowGroupRemoval <- FALSE
probSeqsPrev <- integer() # suspected problem sequences from prior iteration
for (i in seq_len(maxIterations)) {
  cat("Training iteration: ", i, "\n", sep="")
  # train the classifier
  trainingSet <- LearnTaxa(seqs[!remove],
                           names(seqs)[!remove])  # not bothering to make the taxid file, i'm still unsure what it even does (and it is optional)
  # look for problem sequences
  probSeqs <- trainingSet$problemSequences$Index
  if (length(probSeqs)==0) {
    cat("No problem sequences remaining.\n")
    break
  } else if (length(probSeqs)==length(probSeqsPrev) &&
             all(probSeqsPrev==probSeqs)) {
    cat("Iterations converged.\n")
    break
  }
  if (i==maxIterations)
    break
  probSeqsPrev <- probSeqs
  # remove any problem sequences
  index <- which(!remove)[probSeqs]
  remove[index] <- TRUE # remove all problem sequences
  if (!allowGroupRemoval) {
    # replace any removed groups
    missing <- !(u_groups %in% groups[!remove])
    missing <- u_groups[missing]
    if (length(missing) > 0) {
      index <- index[groups[index] %in% missing]
      remove[index] <- FALSE # don't remove
    }
  }
  sum(remove) # total number of sequences eliminated
}

message("A total of ", length(probSeqs), " remain problem seqs after 10 rounds of iterative training")

message("Saving Genus level IDTaxa database")
trainingSet %>% write_rds(paste0(opts$o,"/IdTaxa_20231215.silva.seed_v138_1.ng_FL.RData"))  # again, make sure this is right




