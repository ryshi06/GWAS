

rm(list = ls())
gc()

###--------------------------------------------------------------------------###
### Find Pairs more than 1 Dups
### 
### 20240516
###--------------------------------------------------------------------------###

batches <- c(paste0("GWAS",1:5),"GEM")

for (batch in batches){
  
  # read in the Removal/Within-Batch Dups Info
  Removal <- openxlsx::read.xlsx("/Users/shiruyu/Desktop/Harmonization_Sample_Stats.xlsx", sheet = batch)
  head(Removal)
  if (nrow(Removal[Removal$Removal.Reasons == "Within-batch Duplicates by genotype",]) > 0){
    Dups <- openxlsx::read.xlsx("/Users/shiruyu/Desktop/Harmonization_Sample_Stats.xlsx", sheet = paste0(batch," within batch dups"))
    # head(Dups)
  } else {
    print(paste0("No removal found for batch ", batch, "."))
    next  # skip to the next iteration of the loop
  }
  
  # generate a new data frame containing the Dup Pairs Info
  Dup_Pair <- data.frame(ID1 = intersect(Removal[Removal$Removal.Reasons == "Within-batch Duplicates by genotype",]$SampleID,
                                         unique(c(Dups$Ind1, Dups$Ind2))),
                         ID2 = NA, 
                         ID3 = NA)
  
  for (ind in Dup_Pair$ID1){
    
    if (length(grep(ind, Dups$Ind1)) > 0){
      tmp <- Dups[Dups$Ind1 == ind,]
      IDs <- unique(tmp$Ind2)
    } else {
      tmp <- Dups[Dups$Ind2 == ind,]
      IDs <- unique(tmp$Ind1)
    }
    
    if (length(IDs) > 1){
      Dup_Pair[Dup_Pair$ID1 == ind, ]$ID2 <- IDs[1]
      Dup_Pair[Dup_Pair$ID1 == ind, ]$ID3 <- IDs[2]
    } else {
      Dup_Pair[Dup_Pair$ID1 == ind, ]$ID2 <- IDs[1]
      Dup_Pair[Dup_Pair$ID1 == ind, ]$ID3 <- NA
    }
    
  }
  
  # check if the IDs match between two lists
  Removal_WithinBatch <- sort(Removal[Removal$Removal.Reasons == "Within-batch Duplicates by genotype",]$SampleID)
  Dup_Ind1 <- sort(Dup_Pair$ID1)
 
  if (isTRUE(all.equal(Removal_WithinBatch, Dup_Ind1))){
    
    print(paste0(length(Removal_WithinBatch)," pairs of within-batch duplicates found for batch ", batch, "."))
    # merge the pairs together
    Dup_Pair$Pairs <- ifelse(is.na(Dup_Pair$ID3), paste0(Dup_Pair$ID1, "/", Dup_Pair$ID2), 
                             paste0(Dup_Pair$ID1, "/", Dup_Pair$ID2, "/", Dup_Pair$ID3))
    
    # View(Dup_Pair)
    updated_Removal <- merge(Removal, Dup_Pair[,c("ID1","Pairs")], 
                                    by.x = "SampleID", by.y = "ID1", all.x = TRUE) %>%
      arrange(Removal.Reasons)
    
    # update the final excel 
    
    library(XLConnect)
    writeWorksheetToFile(file = "/Users/shiruyu/Desktop/Harmonization_Sample_Stats_20240516.xlsx", data = updated_Removal, sheet = batch)
    
  } else {
    stop("Error: Removal_WithinBatch and Dup_Ind1 do not match.")
  }
  
}

# [1] "No removal found for batch GWAS1."
# [1] "7 pairs of within-batch duplicates found for batch GWAS2."
# [1] "2 pairs of within-batch duplicates found for batch GWAS3."
# [1] "85 pairs of within-batch duplicates found for batch GWAS4."
# [1] "3 pairs of within-batch duplicates found for batch GWAS5."
# [1] "51 pairs of within-batch duplicates found for batch GEM."
