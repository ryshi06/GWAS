
rm(list = ls())
gc()

library(data.table)
library(dplyr)

###--------------------------------------------------------------------------###
### load Plink IBD
###--------------------------------------------------------------------------###

ibd <- fread("/ix/kfan/Ruyu/harmonization/All_Cohorts/QC_Check/AllCohorts_AllChr_IBD.genome")
nrow(ibd[ibd$PI_HAT > 0.9,])
# 148

ibd_threshold <- ibd[ibd$PI_HAT > 0.9,]

###--------------------------------------------------------------------------###
### load .imiss
###--------------------------------------------------------------------------###

imiss <- fread("/ix/kfan/Ruyu/harmonization/All_Cohorts/QC_Check/AllCohorts_AllChr_Annot_Missing.imiss")
# all equals to 0

###--------------------------------------------------------------------------###
### update Batch info
###--------------------------------------------------------------------------###

pheno <- fread("/ix/kfan/Ruyu/harmonization/All_Cohorts/QC_Check/All_Cohorts_2nd_Imputation_Phenotype.txt")

ibd1 <- merge(ibd_threshold, pheno[,c("scanName","Batch")], by.x = "IID1", by.y = "scanName", all.x = TRUE)
colnames(ibd1) <- c("IID1","FID1","FID2","IID2","RT","EZ","Z0","Z1","Z2","PI_HAT","PHE","DST","PPC","RATIO","Batch1")

ibd2 <- merge(ibd1, pheno[,c("scanName","Batch")], by.x = "IID2", by.y = "scanName", all.x = TRUE)
colnames(ibd2) <- c("IID2","IID1","FID1","FID2","RT","EZ","Z0","Z1","Z2","PI_HAT","PHE","DST","PPC","RATIO","Batch1","Batch2")

# set priority order
priority_order <- c("GWAS5", "GEM", "GWAS4", "GWAS3", "GWAS2", "GWAS1")

determine_keep <- function(val1, val2, priority) {
  # find the positions of val1 and val2 in the priority order
  pos1 <- match(val1, priority)
  pos2 <- match(val2, priority)
  
  # return the value with the higher priority (lower position number)
  if (pos1 < pos2) {
    return(val1)
  } else {
    return(val2)
  }
}

ibd3 <- ibd2 %>%
  rowwise() %>%
  mutate(Keep = determine_keep(Batch1, Batch2, priority_order)) %>%
  ungroup()

table(ibd3$Batch1, ibd3$Keep)
table(ibd3$Batch2, ibd3$Keep)

determine_remove <- function(keep, val1, val2, id1, id2) {
  if (keep == val1) {
    return(id2)
  } else if (keep == val2) {
    return(id1)
  } else {
    return(NA)
  }
}

# apply the function to create the "Remove" column
ibd3 <- ibd3 %>%
  rowwise() %>%
  mutate(Remove = determine_remove(Keep, Batch1, Batch2, IID1, IID2)) %>%
  ungroup()

ibd4 <- ibd3[,c("IID1","Batch1","IID2","Batch2","Keep","Remove")]
ibd4$RemoveBatch <- ifelse(ibd4$Remove == ibd4$IID1, ibd4$Batch1, 
                           ifelse(ibd4$Remove == ibd4$IID2, ibd4$Batch2, NA))

# check duplicates
length(ibd4$Remove[duplicated(ibd4$Remove)])

dups <- ibd4$Remove[duplicated(ibd4$Remove)]
dupIDs <- NA
for(dup in dups){
  tmp <- ibd4[ibd4$Remove == dup,]
  dupIDs <- c(dupIDs, unique(tmp$IID1), unique(tmp$IID2))
}

dupIDs <- unique(dupIDs)[!is.na(unique(dupIDs))]
ibd5 <- rbind(ibd4[ibd4$IID1 %in% dupIDs,],
              ibd4[ibd4$IID2 %in% dupIDs,])

ibd_removal <- ibd3[,c("Remove","Remove")] %>% distinct()

write.table(ibd_removal, "/Users/shiruyu/Desktop/IBD_Remove_list.txt", 
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

table(ibd4$Batch1, ibd4$Keep)

###--------------------------------------------------------------------------###
### check rest Batch info
###--------------------------------------------------------------------------###

pheno1 <- pheno[!(pheno$scanName %in% ibd_removal$Remove), ]
table(pheno1$Batch)

#  GEM GWAS1 GWAS2 GWAS3 GWAS4 GWAS5 
# 2731  2275   619   493  4258   235 

###--------------------------------------------------------------------------###
### check IBD again
###--------------------------------------------------------------------------###

ibd2 <- fread("/ix/kfan/Ruyu/harmonization/All_Cohorts/QC_Check/AllCohorts_AllChr_IBD2.genome")
nrow(ibd2[ibd2$PI_HAT > 0.9,])
