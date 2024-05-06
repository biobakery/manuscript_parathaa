

performance.table <- function(compareData, level){
  
  #If we are comparing species level data
  if(level=="Species"){
    
    #We treat true positives as anything that we set the FLAG variable as true
    TP.parathaa <- compareData %>% filter(Flag.y) %>% nrow()
    #We consider false positives to be anything that the FLAG variable is set to false
    FP.parathaa <- compareData %>% filter(!Flag.y) %>% nrow()
    #We do not consider true negatives in our analysis
    TN.parathaa <- 0
    #False negatives is cases where taxonomy is not assigned and the flag variable is not set
    FN.parathaa <- compareData %>% filter(is.na(Species.parathaa) & is.na(Flag.y)) %>% nrow()
    
    #multi correct are cases where a species is set but it doesn't equal the silva species and the flag is true
    multCorrect.parathaa <- compareData %>% filter(!Species.parathaa==Species.silva & Flag.y) %>% nrow() /nrow(compareData)
    #unique correct are cases when taxonomy match up
    uniqueCorrect.parathaa <- compareData %>% filter(Species.parathaa==Species.silva) %>% nrow() /nrow(compareData)
    #unassigned correct are cases where we think not assigning a taxonomy is the correct choice
    #this occurs when the reference query is from a taxonomy that is outside of the database's knowledge
    unassignedCorrect.parathaa <- compareData %>% filter(is.na(Species.parathaa) & Flag.y) %>% nrow() /nrow(compareData)
    #Unassigned incorrect are cases when we don't make an assignment but the taxonomy we were suppose to assign does exist 
    # within the database
    unassignedIncorrect.parathaa <- compareData %>% filter(is.na(Species.parathaa) & is.na(Flag.y)) %>% nrow() /nrow(compareData)
  }
  
  # For information about individual variables see above section
  if(level=="Genus"){
    TP.parathaa <- compareData %>% filter( Flag.genus.y) %>% nrow()
    FP.parathaa <- compareData %>% filter(!Flag.genus.y) %>% nrow()
    TN.parathaa <- 0
    FN.parathaa <- compareData %>% filter(is.na(Genus.parathaa) & is.na(Flag.genus.y)) %>% nrow()
    
    multCorrect.parathaa <- compareData %>% filter(!Genus.parathaa==Genus.silva & Flag.genus.y) %>% nrow() /nrow(compareData)
    uniqueCorrect.parathaa <- compareData %>% filter(Genus.parathaa==Genus.silva) %>% nrow() /nrow(compareData)
    unassignedCorrect.parathaa <- compareData %>% filter(is.na(Genus.parathaa) & Flag.genus.y) %>% nrow() /nrow(compareData)
    unassignedIncorrect.parathaa <- compareData %>% filter(is.na(Genus.parathaa) & is.na(Flag.genus.y)) %>% nrow() /nrow(compareData)
  }
  
  #Calculate metrics
  accuracy.parathaa <- (TP.parathaa + TN.parathaa) / (TP.parathaa + TN.parathaa + FP.parathaa + FN.parathaa)
  precision.parathaa <- TP.parathaa / (TP.parathaa + FP.parathaa)
  recall.parathaa <- TP.parathaa / (TP.parathaa + FN.parathaa)
  f1.parathaa <- 2 * (precision.parathaa * recall.parathaa) / (precision.parathaa + recall.parathaa)
  fpr.parathaa <- FP.parathaa / nrow(compareData)
  
  #Do the same as above but for dada2
  if(level=="Species"){
    TP.dada <- compareData %>% filter( Flag.x) %>% nrow()
    FP.dada <- compareData %>% filter(!Flag.x) %>% nrow()
    TN.dada <- 0
    FN.dada <- compareData %>% filter(is.na(Species.dada) & is.na(Flag.x)) %>% nrow()
    multCorrect.dada <- compareData %>% filter(!Species.dada==Species.silva & Flag.x) %>% nrow() /nrow(compareData)
    uniqueCorrect.dada <- compareData %>% filter(Species.dada==Species.silva) %>% nrow() /nrow(compareData)
    unassignedCorrect.dada <- compareData %>% filter(is.na(Species.dada) & Flag.x) %>% nrow() /nrow(compareData)
    unassignedIncorrect.dada <- compareData %>% filter(is.na(Species.dada) & is.na(Flag.x)) %>% nrow() /nrow(compareData)
  }
  if(level=="Genus"){
    TP.dada <- compareData %>% filter( Flag.genus.x) %>% nrow()
    FP.dada <- compareData %>% filter(!Flag.genus.x) %>% nrow()
    TN.dada <- 0
    FN.dada <- compareData %>% filter(is.na(Genus.dada) & is.na(Flag.genus.x)) %>% nrow()
    multCorrect.dada <- compareData %>% filter(!Genus.dada==Genus.silva & Flag.genus.x) %>% nrow() /nrow(compareData)
    uniqueCorrect.dada <- compareData %>% filter(Genus.dada==Genus.silva) %>% nrow() /nrow(compareData)
    unassignedCorrect.dada <- compareData %>% filter(is.na(Genus.dada) & Flag.genus.x) %>% nrow() /nrow(compareData)
    unassignedIncorrect.dada <- compareData %>% filter(is.na(Genus.dada) & is.na(Flag.genus.x)) %>% nrow() /nrow(compareData)
  }
  
  accuracy.dada <- (TP.dada + TN.dada) / (TP.dada + TN.dada + FP.dada + FN.dada)
  precision.dada <- TP.dada / (TP.dada + FP.dada)
  recall.dada <- TP.dada / (TP.dada + FN.dada)
  f1.dada <- 2 * (precision.dada * recall.dada) / (precision.dada + recall.dada)
  fpr.dada <- FP.dada / nrow(compareData)

  
  #set up the output table
  rows1 <- c(             "Accuracy", "Precision", "Recall", "F1 Score",
                          "Uniquely Correct", "One-to-many Correct", "Incorrect", "Unassigned Correct", "Unassigned Incorrect")
  
  table.out <- matrix(NA, nrow = length(rows1), ncol=2)
  colnames(table.out) <- c("Parathaa", "DADA2")
  rownames(table.out) <- rows1
  table.out["Uniquely Correct", "Parathaa"] <- uniqueCorrect.parathaa
  table.out["Uniquely Correct", "DADA2"] <- uniqueCorrect.dada
  table.out["One-to-many Correct", "Parathaa"] <- multCorrect.parathaa
  table.out["One-to-many Correct", "DADA2"] <- multCorrect.dada
  table.out["Incorrect", "Parathaa"] <- fpr.parathaa
  table.out["Incorrect", "DADA2"] <- fpr.dada
  table.out["Unassigned Correct", "Parathaa"] <- unassignedCorrect.parathaa
  table.out["Unassigned Correct", "DADA2"] <- unassignedCorrect.dada
  table.out["Unassigned Incorrect", "Parathaa"] <- unassignedIncorrect.parathaa
  table.out["Unassigned Incorrect", "DADA2"] <- unassignedIncorrect.dada
  
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
