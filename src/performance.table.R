

performance.table <- function(compareData, level){
  
  if(level=="Species"){
    TP.parathaa <- compareData %>% filter(Flag.y) %>% nrow()
    FP.parathaa <- compareData %>% filter(!Flag.y) %>% nrow()
    TN.parathaa <- 0
    FN.parathaa <- compareData %>% filter(is.na(Species.parathaa)) %>% nrow()
    
    multCorrect.parathaa <- compareData %>% filter(!Species.parathaa==Species.silva & Flag.y) %>% nrow() /nrow(compareData)
    uniqueCorrect.parathaa <- compareData %>% filter(Species.parathaa==Species.silva) %>% nrow() /nrow(compareData)
    unassigned.parathaa <- compareData %>% filter(is.na(Flag.y)) %>% nrow() /nrow(compareData)
  }
  
  if(level=="Genus"){
    TP.parathaa <- compareData %>% filter( Flag.genus.y) %>% nrow()
    FP.parathaa <- compareData %>% filter(!Flag.genus.y) %>% nrow()
    TN.parathaa <- 0
    FN.parathaa <- compareData %>% filter(is.na(Genus.parathaa)) %>% nrow()
    
    multCorrect.parathaa <- compareData %>% filter(!Genus.parathaa==Genus.silva & Flag.genus.y) %>% nrow() /nrow(compareData)
    uniqueCorrect.parathaa <- compareData %>% filter(Genus.parathaa==Genus.silva) %>% nrow() /nrow(compareData)
    unassigned.parathaa <- compareData %>% filter(is.na(Flag.genus.y)) %>% nrow() /nrow(compareData)
  }
  
  accuracy.parathaa <- (TP.parathaa + TN.parathaa) / (TP.parathaa + TN.parathaa + FP.parathaa + FN.parathaa)
  precision.parathaa <- TP.parathaa / (TP.parathaa + FP.parathaa)
  recall.parathaa <- TP.parathaa / (TP.parathaa + FN.parathaa)
  f1.parathaa <- 2 * (precision.parathaa * recall.parathaa) / (precision.parathaa + recall.parathaa)
  fpr.parathaa <- FP.parathaa / nrow(compareData)

  if(level=="Species"){
    TP.dada <- compareData %>% filter( Flag.x) %>% nrow()
    FP.dada <- compareData %>% filter(!Flag.x) %>% nrow()
    TN.dada <- 0
    FN.dada <- compareData %>% filter(is.na(Species.dada)) %>% nrow()
    multCorrect.dada <- compareData %>% filter(!Species.dada==Species.silva & Flag.x) %>% nrow() /nrow(compareData)
    uniqueCorrect.dada <- compareData %>% filter(Species.dada==Species.silva) %>% nrow() /nrow(compareData)
    unassigned.dada <- compareData %>% filter(is.na(Flag.x)) %>% nrow() /nrow(compareData)
  }
  if(level=="Genus"){
    TP.dada <- compareData %>% filter( Flag.genus.x) %>% nrow()
    FP.dada <- compareData %>% filter(!Flag.genus.x) %>% nrow()
    TN.dada <- 0
    FN.dada <- compareData %>% filter(is.na(Genus.dada)) %>% nrow()
    multCorrect.dada <- compareData %>% filter(!Genus.dada==Genus.silva & Flag.genus.x) %>% nrow() /nrow(compareData)
    uniqueCorrect.dada <- compareData %>% filter(Genus.dada==Genus.silva) %>% nrow() /nrow(compareData)
    unassigned.dada <- compareData %>% filter(is.na(Flag.genus.x)) %>% nrow() /nrow(compareData)
  }
  
  accuracy.dada <- (TP.dada + TN.dada) / (TP.dada + TN.dada + FP.dada + FN.dada)
  precision.dada <- TP.dada / (TP.dada + FP.dada)
  recall.dada <- TP.dada / (TP.dada + FN.dada)
  f1.dada <- 2 * (precision.dada * recall.dada) / (precision.dada + recall.dada)
  fpr.dada <- FP.dada / nrow(compareData)

  
  rows1 <- c(             "Accuracy", "Precision", "Recall", "F1 Score",
                          "Uniquely Correct", "One-to-many Correct", "Incorrect", "Unassigned")
  
  table.out <- matrix(NA, nrow = length(rows1), ncol=2)
  colnames(table.out) <- c("Parathaa", "DADA2")
  rownames(table.out) <- rows1
  table.out["Uniquely Correct", "Parathaa"] <- uniqueCorrect.parathaa
  table.out["Uniquely Correct", "DADA2"] <- uniqueCorrect.dada
  table.out["One-to-many Correct", "Parathaa"] <- multCorrect.parathaa
  table.out["One-to-many Correct", "DADA2"] <- multCorrect.dada
  table.out["Incorrect", "Parathaa"] <- fpr.parathaa
  table.out["Incorrect", "DADA2"] <- fpr.dada
  table.out["Unassigned", "Parathaa"] <- unassigned.parathaa
  table.out["Unassigned", "DADA2"] <- unassigned.dada
  
  table.out["Accuracy", "Parathaa"] <- accuracy.parathaa
  table.out["Accuracy", "DADA2"] <- accuracy.dada
  table.out["Precision", "Parathaa"] <- precision.parathaa
  table.out["Precision", "DADA2"] <- precision.dada
  table.out["Recall", "Parathaa"] <- recall.parathaa
  table.out["Recall", "DADA2"] <- recall.dada
  table.out["F1 Score", "Parathaa"] <- f1.parathaa
  table.out["F1 Score", "DADA2"] <- f1.dada


  print(round(table.out, 3))
}

if(FALSE){
### Shortcut to get tables without having to rerun Plots.Table.1.R
load("output/Figures/test_2e12/V1V2_full_comparisons.RData")
outputDir <- "output/Figures/test_2e12/"
regionName <- "V1V2"
t1 <- performance.table(compare.synth, "Species")
write.table(t1, file= file.path(outputDir, paste0(regionName, "_Species_performance.tsv")), sep="\t", col.names = NA)
t2 <- performance.table(compare.synth, "Genus")
write.table(t2, file= file.path(outputDir, paste0(regionName, "_Genus_performance.tsv")), sep="\t", col.names = NA)


load("output/Figures/test_2e12/V4V5_full_comparisons.RData")
outputDir <- "output/Figures/test_2e12/"
regionName <- "V4V5"
t1 <- performance.table(compare.synth, "Species")
write.table(t1, file= file.path(outputDir, paste0(regionName, "_Species_performance.tsv")), sep="\t", col.names = NA)
t2 <- performance.table(compare.synth, "Genus")
write.table(t2, file= file.path(outputDir, paste0(regionName, "_Genus_performance.tsv")), sep="\t", col.names = NA)

load("output/Figures/synth_mult/V1V2_full_comparisons.RData")
performance.table(compare.synth, "Species")

load("output/Figures/synth_mult/V4V5_full_comparisons.RData")
performance.table(compare.synth, "Species")

}
