
###--------------------------------------------------------------------------###
### clean up the Sex and Race in Harmonization Samples
###--------------------------------------------------------------------------###

rm(list = ls())
setwd("/ix/kfan/Ruyu/harmonization/All_Cohorts/Ref_files_2024_06_14")
# update reference file's sex and race

library(dplyr)
library(openxlsx)

### load the old reference 
oldRef <- openxlsx::read.xlsx("V_DataRequest_Ruyu_GWAS_Info_2023_12_20.xlsx")
oldRef$Scan_Name2 <- gsub(" ", "", oldRef$Scan_Name)

oldRef$sex_recode <- ifelse(!is.na(oldRef$Sex_DB) & oldRef$Sex_DB > 0, oldRef$Sex_DB, oldRef$Sex)
oldRef$sex_recode1 <- ifelse(oldRef$sex_recode < 0, NA, oldRef$sex_recode)

# table(oldRef$Race)
# sum(is.na(oldRef$Race))
oldRef$race_recode <- ifelse(oldRef$Race %in% c("White", "W"), "White", 
                             ifelse(oldRef$Race %in% c("B", "Black"), "Black", 
                                    ifelse(oldRef$Race %in% c("A"), "Asian",
                                           ifelse(oldRef$Race < 0, NA, "Others"))))
table(oldRef$race_recode)
sum(is.na(oldRef$race_recode))

# table(oldRef$Race_DB)
# sum(is.na(oldRef$Race_DB))

oldRef$race_db_recode <- ifelse(oldRef$Race_DB %in% c("White", "W"), "White", 
                                ifelse(oldRef$Race_DB %in% c("B", "Black"), "Black", 
                                       ifelse(oldRef$Race_DB %in% c("A"), "Asian",
                                              ifelse(oldRef$Race_DB < 0, NA, "Others"))))
table(oldRef$race_db_recode)
sum(is.na(oldRef$race_db_recode))

oldRef$race_recode1 <- ifelse(!is.na(oldRef$race_db_recode), oldRef$race_db_recode, oldRef$race_recode)
table(oldRef$race_recode1)

oldRef$V1 <- paste0(oldRef$scanName, "_", oldRef$scanName)

### load GWAS4_clinical_updated_info
GWAS4_Info <- read.csv("GWAS_4_Clinical_Info_updated_2023_12_13.csv")

GWAS4_Info$V2 <- gsub(" ", "", GWAS4_Info$V1)
table(GWAS4_Info$Race)

GWAS4_Info$race_recode <- ifelse(GWAS4_Info$Race %in% c("White", "W"), "White", 
                                 ifelse(GWAS4_Info$Race %in% c("B", "Black"), "Black", 
                                        ifelse(GWAS4_Info$Race %in% c("A"), "Asian",
                                               ifelse(GWAS4_Info$Race < 0, NA, "Others"))))

table(GWAS4_Info$race_recode)
sum(is.na(GWAS4_Info$race_recode))

length(intersect(GWAS4_Info$V1, oldRef$Scan_Name))
# [1] 4539

### replace the gwas4 part using the updated info

GWAS_Ref <- oldRef[!(oldRef$Scan_Name %in% intersect(GWAS4_Info$V1, oldRef$Scan_Name)), 
                   c("AlzID","Scan_Name","LabID","StudyID","Study","sex_recode1","race_recode1")] %>% distinct()
dim(GWAS_Ref)
colnames(GWAS_Ref) <- c("AlzID","V1","LabID","StudyID","Study","Sex","Race")

updatedRef <- rbind(GWAS_Ref[,c("AlzID","V1","Study","Sex","Race")], 
                    GWAS4_Info[,c("AlzID","V1","Study","Sex","Race")])

setdiff(updatedRef$V1,oldRef$Scan_Name)
setdiff(oldRef$Scan_Name,updatedRef$V1)

length(unique(updatedRef$V1))

# recode sex and race
updatedRef$sex_recode <- ifelse(updatedRef$Sex < 0, NA, updatedRef$Sex)
table(updatedRef$sex_recode)
#    F    M 
# 6283 4969
sum(is.na(updatedRef$sex_recode))
# 80
updatedRef[is.na(updatedRef$sex_recode), "sex_recode"] <- "F"

updatedRef$race_recode <- ifelse(updatedRef$Race %in% c("White", "W"), "White", 
                                 ifelse(updatedRef$Race %in% c("B", "Black"), "Black", 
                                        ifelse(updatedRef$Race %in% c("A", "Asian"), "Asian",
                                               ifelse(updatedRef$Race < 0, NA, "Others"))))
table(updatedRef$race_recode)
# Asian  Black Others  White 
#    27    707     74   9117
sum(is.na(updatedRef$race_recode))
# 1407

## write.xlsx(updatedRef, "V_DataRequest_Ruyu_GWAS_Info_2023_12_20_updatedGWAS4.xlsx")

### 2nd Phenotype
pheno <- fread("../All_Cohorts_2nd_Imputation_Phenotype.txt")

updatedRef$V2 <- gsub(" ", "", updatedRef$V1)
for (seq in 1:nrow(updatedRef)){
  if(updatedRef[seq,]$V2 %in% pheno$scanName){
    updatedRef[seq, "sex_recode"] <- unique(pheno[pheno$scanName %in% updatedRef[seq,]$V2, ]$sex_recode)
  }
}

### load all samples N = 11387
### batch specific sample name
batches1 <- c(paste0("GWAS",1:5), "GEM")
batches2 <- c(paste0("GWAS_",1:5), "GWAS_GEM")

ori_list <- list()
for (batch in batches2){
  ori_list[[batch]] <- read.csv(paste0(batch, "_2023_07_24.csv"))
  ori_list[[batch]]$Batch <- batch
}
ori <- do.call(rbind, ori_list)
# 11113

noDups_list <- list()
for (batch in batches1){
  noDups_list[[batch]] <- fread(paste0(batch, "_plink_IDs_noDup.txt"))
  noDups_list[[batch]]$Batch <- batch
}
noDups <- do.call(rbind, noDups_list)
# 10816

allScan_list <- list()
for (batch in batches1){
  allScan_list[[batch]] <- fread(paste0(batch, "_allChr_scanName.txt"), header = FALSE)
  allScan_list[[batch]]$Batch <- batch
}
allScan <- do.call(rbind, allScan_list)
# 11020
allScan$V1 <- substr(allScan$V1, 1, nchar(allScan$V1)/2)

ex_list <- list()
for (batch in batches1){
  ex_list[[batch]] <- openxlsx::read.xlsx("Harmonization_Sample_Stats.xlsx", sheet = batch)
  if(batch != "GWAS1"){
    ex_list[[batch]]$Batch <- batch
  }
}
ex <- do.call(rbind, ex_list)
rownames(ex) <- NULL
unique(ex$Removal.Reasons)
ex1 <- ex[ex$Removal.Reasons %in% c("Missing call rate > 5%", "Within-batch Duplicates by genotype"),]
colnames(ex1) <- c("V1","Removal.Reasons","Pairs","Batch")
sex.ex <- ex[ex$Removal.Reasons == "mislabeled sex",]
cross.batch.dups.ex <- ex[ex$Removal.Reasons == "Cross-batch Duplicates by overlapped Sample ID across batches",]

gwas1.fam <- fread("Omni1-quad_07_01_2010_GenTrain_2_Kamboh_case_control.fam")
gwas1.fam$scanName <- paste0(gwas1.fam$V1, "_", gwas1.fam$V2)
gwas1.fam$scanName2 <- paste0(gwas1.fam$scanName, "_", gwas1.fam$scanName)

gwas1.ex <- as.data.frame(setdiff(gwas1.fam$scanName, allScan[allScan$Batch == "GWAS1",]$V1))
head(gwas1.ex)
gwas1.ex$`Removal.Reasons` <- "Missing call rate > 5%"
colnames(gwas1.ex) <- c("SampleID", "Removal.Reasons")
gwas1.ex$Batch <- "GWAS1"
# gwas1.ex$V1 <- substr(gwas1.ex$SampleID, 1, nchar(gwas1.ex$SampleID)/2)
gwas1.ex$V2 <- sub("^[^_]*_", "", gwas1.ex$SampleID)

colnames(gwas1.ex) <- c("V1", "Removal.Reasons", "Batch", "V2")
ex2 <- rbind(gwas1.ex[,c("V1", "Removal.Reasons", "Batch")], ex1[,c("V1", "Removal.Reasons", "Batch")])
dim(ex2)

###--------------------------------------------------------------------------###
### get the total N = 11387 list with GWAS QC Pass/Fail info
###--------------------------------------------------------------------------###

allScan$Removal.Reasons <- as.character(NA_character_)
all <- rbind(ex2, allScan)
table(all$Batch)
#  GEM GWAS1 GWAS2 GWAS3 GWAS4 GWAS5 
# 2800  2440   720   640  4539   248
# View(all)

# fill in sex mislabeled info
sex_ex_samples <- intersect(sex.ex$SampleID, all$V1)

for (sample in sex_ex_samples) {
  
  sex.ex_tmp <- sex.ex[sex.ex$SampleID %in% sample,]
  sex.ex_tmp.batch <- unique(sex.ex_tmp$Batch)
  
  if (is.na(all[all$V1 %in% sample & all$Batch %in% sex.ex_tmp.batch, ]$Removal.Reasons)){
    all[all$V1 %in% sample & all$Batch %in% sex.ex_tmp.batch, "Removal.Reasons"] <- "Mislabeled Sex"
  }
  
}

# fill in cross-batch sample ID overlap info
dups_ex_samples <- intersect(cross.batch.dups.ex$SampleID, all$V1)

for (sample in dups_ex_samples) {
  dups_tmp <- cross.batch.dups.ex[cross.batch.dups.ex$SampleID %in% sample,]
  dups_tmp.batch <- unique(dups_tmp$Batch)
  
  if(nrow(dups_tmp) > 1){
    batch1 <- unique(dups_tmp$Batch)[1]
    batch2 <- unique(dups_tmp$Batch)[2]
    
    if (is.na(all[all$V1 %in% sample & all$Batch %in% batch1, ]$Removal.Reasons)){
      all[all$V1 %in% sample & all$Batch %in% batch1, "Removal.Reasons"] <- "Cross-batch Duplicates by overlapped Sample ID across batches"
    }
    
    if (is.na(all[all$V1 %in% sample & all$Batch %in% batch2, ]$Removal.Reasons)){
      all[all$V1 %in% sample & all$Batch %in% batch2, "Removal.Reasons"] <- "Cross-batch Duplicates by overlapped Sample ID across batches"
    }
    
  } else {
    if (is.na(all[all$V1 %in% sample & all$Batch %in% dups_tmp.batch, ]$Removal.Reasons)){
      all[all$V1 %in% sample & all$Batch %in% dups_tmp.batch, "Removal.Reasons"] <- "Cross-batch Duplicates by overlapped Sample ID across batches"
    }
  }
}

### assign GWAS QC Pass/Fail to Samples
all$GWAS_QC <- ifelse(is.na(all$Removal.Reasons), "Pass", "Failed")
table(all$GWAS_QC)
# Failed   Pass 
# 634  10753

table(all[is.na(all$GWAS_QC),]$Batch)
# sum(is.na(all$Removal.Reasons))
table(all$Removal.Reasons)

### get the final Ref file
all$Sex <- as.character(NA_character_)
all$Race <- as.character(NA_character_)

for (seq in 1:nrow(all)){
  if(all[seq,]$V1 %in% updatedRef$V2){
    all[seq, "Sex"] <- unique(updatedRef[updatedRef$V2 %in% all[seq,]$V1, ]$sex_recode)
    all[seq, "Race"] <- unique(updatedRef[updatedRef$V2 %in% all[seq,]$V1, ]$race_recode)
  }
}

# View(all[is.na(all$Sex),])
nrow(all[is.na(all$Sex),])

# length(intersect(all[is.na(all$Sex),]$V1, oldRef$scanName))

for (seq in 1:nrow(all)) {
  if (is.na(all$Sex[seq])) {
    matching_row <- pheno[pheno$scanName == all$V1[seq], ]
    if (nrow(matching_row) > 0) {
      all$Sex[seq] <- unique(matching_row$sex_recode)
    }
  }
}

oldRef$V1 <- paste0(oldRef$Scan_Name, "_", oldRef$Scan_Name)

for (seq in 1:nrow(all)) {
  if (is.na(all$Sex[seq])) {
    matching_row <- oldRef[oldRef$V1 == all$V1[seq], ]
    if (nrow(matching_row) > 0) {
      all$Sex[seq] <- unique(matching_row$sex_recode)
    }
  }
}

all[is.na(all$Sex), "Sex"] <- "F"
table(all$Sex)
#    F    M 
# 6389 4998 

for (seq in 1:nrow(all)) {
  if (is.na(all$Race[seq])) {
    matching_row <- pheno[pheno$scanName == all$V1[seq], ]
    if (nrow(matching_row) > 0) {
      all$Race[seq] <- unique(matching_row$race_recode)
    }
  }
}

for (seq in 1:nrow(all)) {
  if (is.na(all$Race[seq])) {
    matching_row <- oldRef[oldRef$V1 == all$V1[seq], ]
    if (nrow(matching_row) > 0) {
      all$Race[seq] <- unique(matching_row$race_recode)
    }
  }
}

for (seq in 1:nrow(all)) {
  if (is.na(all$Race[seq])) {
    matching_row <- oldRef[oldRef$Scan_Name == all$V1[seq], ]
    if (nrow(matching_row) > 0) {
      all$Race[seq] <- unique(matching_row$race_recode)
    }
  }
}

table(all$Race)
# View(all[is.na(all$Race),])
all[is.na(all$Race), "Race"] <- "Unknown"

### add AlzID
all$AlzID <- as.character(NA_character_)

# cross check with the APOE24 table
### the APOE24 table
apoe24 <- openxlsx::read.xlsx("APOE24_Sample_Characteristics_SENIORS-updated_N15186_2024_06_11.xlsx")
apoe24_split <- apoe24 %>%
  tidyr::separate_rows(LabIDList, sep = ",") %>%
  as.data.frame()

for (seq in 1:nrow(all)) {
  if (is.na(all$AlzID[seq])) {
    matching_row <- apoe24_split[toupper(apoe24_split$LabIDList) == toupper(gsub("R","",all$V1[seq])), ]
    if (nrow(matching_row) > 0) {
      all$AlzID[seq] <- unique(matching_row$AlzID)
    }
  }
}

for (seq in 1:nrow(all)){
  if(all[seq,]$V1 %in% pheno$scanName){
    all[seq, "AlzID"] <- unique(pheno[pheno$scanName %in% all[seq,]$V1, ]$AlzID)
  }
}

for (seq in 1:nrow(all)) {
  if (is.na(all$AlzID[seq])) {
    matching_row <- oldRef[oldRef$scanName == all$V1[seq], ]
    if (nrow(matching_row) > 0) {
      all$AlzID[seq] <- unique(matching_row$AlzID_DB)
    }
  }
}

for (seq in 1:nrow(all)) {
  if (is.na(all$AlzID[seq])) {
    matching_row <- oldRef[oldRef$V1 == all$V1[seq], ]
    if (nrow(matching_row) > 0) {
      all$AlzID[seq] <- unique(matching_row$AlzID_DB)
    }
  }
}

for (seq in 1:nrow(all)) {
  if (is.na(all$AlzID[seq])) {
    matching_row <- oldRef[oldRef$Scan_Name2 == all$V1[seq], ]
    if (nrow(matching_row) > 0) {
      all$AlzID[seq] <- unique(matching_row$AlzID_DB)
    }
  }
}

# manually assign by split24 table
all[all$V1 == "G0444_Kamboh_C_0856", "AlzID"] <- "ALZ114179"
all[all$V1 == "M0279_Kamboh_C_0276", "AlzID"] <- "ALZ114768"
all[all$V1 == "M0280_Kamboh_C_0277", "AlzID"] <- "ALZ114769"
all[all$V1 == "S275", "AlzID"] <- "ALZ100914"
all[all$V1 == "S313", "AlzID"] <- "ALZ100923"

# View(all[is.na(all$AlzID),])

###--------------------------------------------------------------------------###
### add study info
###--------------------------------------------------------------------------###

all$Study <- as.character(NA_character_)

for (seq in 1:nrow(all)){
  if(all[seq,]$AlzID %in% apoe24_split$AlzID){
    all[seq, "Study"] <- unique(apoe24_split[apoe24_split$AlzID %in% all[seq,]$AlzID, ]$StudyList)
  }
}

for (seq in 1:nrow(all)) {
  if (is.na(all$Study[seq])) {
    matching_row <- oldRef[oldRef$AlzID_DB %in% all$AlzID[seq], ]
    if (nrow(matching_row) > 0) {
      all$Study[seq] <- unique(matching_row$Study_DB)
    }
  }
}

all[all$AlzID %in% "ALZ103319", "Study"] <- "GEM"
all[is.na(all$Study), "Study"] <- "ADRC"
# View(all[is.na(all$Study),])
sum(is.na(all$Study))
# table(all$Study)
##write.xlsx(all, "../QC_Check/GWAS_Harmonization_N=11387_Demographic_Final.xlsx")



