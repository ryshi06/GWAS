
rm(list = ls())

#####################################
### Check cross-cohort duplicates ###
#####################################

setwd("/zfs1/kfan/Ruyu/harmonization_Sep5/All_Cohorts/Cohort_Check/")
cohort_dups <- openxlsx::read.xlsx("duplicated_harmonized_samples_across_cohort_2024_1_4.xlsx")
cohort_dups <- cohort_dups[cohort_dups$LabID_DB != "-1", ]
dim(cohort_dups)
# [1] 575  13

cohort_dups_GWAS5 <- subset(cohort_dups, GWAS5 == "Y")
dim(cohort_dups_GWAS5)
# [1] 22 13

table(cohort_dups_GWAS5$GWAS4)
# criteria 1: for those genotyped in GWAS5, exclude duplicates in other cohorts
ex1 <- cohort_dups_GWAS5$LabID_DB # all in GWAS4

cohort_dups_oth <- cohort_dups[!(cohort_dups$AlzID_DB %in% cohort_dups_GWAS5$AlzID_DB),]
dim(cohort_dups_oth)
# [1] 553  13

table(cohort_dups_oth$GWAS5)
table(cohort_dups_oth$GWAS4) # 162
table(cohort_dups_oth$GWAS3) # 307
table(cohort_dups_oth$GWAS2) # 119

# integrate missing.e1 information from each cohort
load("/zfs1/kfan/QC_RData/GEM_GWAS_QC.RData")
scanAnnot2$scanName <- toupper(scanAnnot2$scanName)
cohort_dups_missing_e1 <- merge(cohort_dups_oth, scanAnnot2[,c("scanName", "missing.e1")], 
                                by.x = "LabID_DB", by.y = "scanName", all.x = TRUE)
dim(cohort_dups_missing_e1)
# [1] 553  14
colnames(cohort_dups_missing_e1) <- c("LabID_DB","GWAS1","GWAS2","GWAS3","GWAS4","GWAS5","GEM",
                                      "AlzID_DB","StudyID_DB","Sex_DB","Flag_Sex","Race_DB","Flag_Race", 
                                      "missing.e1_GEM")

table(cohort_dups_missing_e1$GWAS5)
table(cohort_dups_missing_e1$GWAS4) # 162
table(cohort_dups_missing_e1$GWAS3) # 307
table(cohort_dups_missing_e1$GWAS2) # 119

load("/zfs1/kfan/QC_RData/GWAS_2_QC.RData")
scanAnnot2$scanName <- toupper(scanAnnot2$scanName)
cohort_dups_missing_e1 <- merge(cohort_dups_missing_e1, scanAnnot2[,c("scanName", "missing.e1")], 
                                by.x = "LabID_DB", by.y = "scanName", all.x = TRUE)
dim(cohort_dups_missing_e1)
# [1] 553  15
colnames(cohort_dups_missing_e1) <- c("LabID_DB","GWAS1","GWAS2","GWAS3","GWAS4","GWAS5","GEM",
                                      "AlzID_DB","StudyID_DB","Sex_DB","Flag_Sex","Race_DB","Flag_Race", 
                                      "missing.e1_GEM", "missing.e1_GWAS2")

table(cohort_dups_missing_e1$GWAS5)
table(cohort_dups_missing_e1$GWAS4) # 162
table(cohort_dups_missing_e1$GWAS3) # 307
table(cohort_dups_missing_e1$GWAS2) # 119

load("/zfs1/kfan/QC_RData/GWAS_3_QC.RData")
scanAnnot2$scanName <- toupper(scanAnnot2$scanName)
cohort_dups_missing_e1 <- merge(cohort_dups_missing_e1, scanAnnot2[,c("scanName", "missing.e1")], 
                                by.x = "LabID_DB", by.y = "scanName", all.x = TRUE)
dim(cohort_dups_missing_e1)
# [1] 525  16
colnames(cohort_dups_missing_e1) <- c("LabID_DB","GWAS1","GWAS2","GWAS3","GWAS4","GWAS5","GEM",
                                      "AlzID_DB","StudyID_DB","Sex_DB","Flag_Sex","Race_DB","Flag_Race", 
                                      "missing.e1_GEM", "missing.e1_GWAS2", "missing.e1_GWAS3")

table(cohort_dups_missing_e1$GWAS5)
table(cohort_dups_missing_e1$GWAS4) # 162
table(cohort_dups_missing_e1$GWAS3) # 307
table(cohort_dups_missing_e1$GWAS2) # 119

load("/zfs1/kfan/Ruyu/harmonization_Sep5/GWAS4/QC_check/GWAS_4_QC.RData")
scanAnnot@data$scanName <- toupper(scanAnnot@data$scanName)
cohort_dups_missing_e1 <- merge(cohort_dups_missing_e1, scanAnnot@data[,c("scanName", "missing.e1")], 
                                by.x = "LabID_DB", by.y = "scanName", all.x = TRUE)
dim(cohort_dups_missing_e1)
# [1] 553  17
colnames(cohort_dups_missing_e1) <- c("LabID_DB","GWAS1","GWAS2","GWAS3","GWAS4","GWAS5","GEM",
                                      "AlzID_DB","StudyID_DB","Sex_DB","Flag_Sex","Race_DB","Flag_Race", 
                                      "missing.e1_GEM", "missing.e1_GWAS2", "missing.e1_GWAS3", "missing.e1_GWAS4")

table(cohort_dups_missing_e1$GWAS5)
table(cohort_dups_missing_e1$GWAS4) # 162
table(cohort_dups_missing_e1$GWAS3) # 307
table(cohort_dups_missing_e1$GWAS2) # 119

# compare values from different columns
library(dplyr)
library(hablar)

tmp <- cohort_dups_missing_e1 %>%
  rowwise() %>%
  mutate(min = min_(c_across(c("missing.e1_GEM", "missing.e1_GWAS2", "missing.e1_GWAS3","missing.e1_GWAS4")))) %>%
  mutate(cohort = Reduce(coalesce, across(c("missing.e1_GEM", "missing.e1_GWAS2", "missing.e1_GWAS3","missing.e1_GWAS4"), ~ ifelse(. == min, cur_column(), NA_character_)))) %>%
  as.data.frame()

tmp$cohort_keep <- gsub("missing.e1_", "", tmp$cohort)
table(tmp$cohort_keep)
# GEM GWAS2 GWAS3 GWAS4 
# 187    51   160    89
sum(is.na(tmp$cohort_keep))
# 66
dim(tmp)
# [1] 553  20

cohort_missing_e1_output <- tmp[!is.na(tmp$cohort_keep), ]
dim(cohort_missing_e1_output)
# [1] 487  20
cohort_missing_e1_output$cohort_count <- rowSums(!is.na(cohort_missing_e1_output[,c("GWAS1","GWAS2","GWAS3","GWAS4","GWAS5","GEM")]))
cohort_missing_e1_output_tmp <- cohort_missing_e1_output[cohort_missing_e1_output$cohort_count >= 2, ]
dim(cohort_missing_e1_output_tmp)
# [1] 406  21

# cohort_missing_e1_problematic <- tmp[is.na(tmp$cohort_keep), ]

ex_gwas1 <- c()
ex_gwas2 <- c()
ex_gwas3 <- c()
ex_gwas4 <- c()
ex_gwas5 <- c()
ex_gem <- c()

for (sample in 1:nrow(cohort_missing_e1_output)){
  
  tmp <- cohort_missing_e1_output[sample, ]
  
  if (!is.na(tmp$GWAS1) & tmp$GWAS1 == "Y" & tmp$cohort != "GWAS1") {
    ex_gwas1 <- c(ex_gwas1, tmp$LabID_DB)
  }
  
  if (!is.na(tmp$GWAS2) & tmp$GWAS2 == "Y" & tmp$cohort != "GWAS2") {
    ex_gwas2 <- c(ex_gwas2, tmp$LabID_DB)
  } 

  if (!is.na(tmp$GWAS3) & tmp$GWAS3 == "Y" & tmp$cohort != "GWAS3") {
    ex_gwas3 <- c(ex_gwas3, tmp$LabID_DB)
  } 
  
  if (!is.na(tmp$GWAS4) & tmp$GWAS4 == "Y" & tmp$cohort != "GWAS4") {
    ex_gwas4 <- c(ex_gwas4, tmp$LabID_DB)
  }
  
  if (!is.na(tmp$GEM) & tmp$GEM == "Y" & tmp$cohort != "GEM") {
    ex_gem <- c(ex_gem, tmp$LabID_DB)
  }
  
}

ex_gwas4 <- c(ex1, ex_gwas4)

length(unique(ex_gwas1))
# 163
length(unique(ex_gwas2))
# 113
length(unique(ex_gwas3))
# 280
length(unique(ex_gwas4))
# 143
length(unique(ex_gem))
# 226

cohort_flag <- cohort_dups[!is.na(cohort_dups$Flag_Race) | !is.na(cohort_dups$Flag_Sex), ]

for (sample in 1:nrow(cohort_flag)){
  
  tmp <- cohort_flag[sample, ]
  
  if (!is.na(tmp$GWAS1) & tmp$GWAS1 == "Y") {
    ex_gwas1 <- c(ex_gwas1, tmp$LabID_DB)
  }
  
  if (!is.na(tmp$GWAS2) & tmp$GWAS2 == "Y") {
    ex_gwas2 <- c(ex_gwas2, tmp$LabID_DB)
  } 
  
  if (!is.na(tmp$GWAS3) & tmp$GWAS3 == "Y") {
    ex_gwas3 <- c(ex_gwas3, tmp$LabID_DB)
  } 
  
  if (!is.na(tmp$GWAS4) & tmp$GWAS4 == "Y") {
    ex_gwas4 <- c(ex_gwas4, tmp$LabID_DB)
  }
  
  if (!is.na(tmp$GWAS5) & tmp$GWAS5 == "Y") {
    ex_gwas5 <- c(ex_gwas5, tmp$LabID_DB)
  }
  
  if (!is.na(tmp$GEM) & tmp$GEM == "Y") {
    ex_gem <- c(ex_gem, tmp$LabID_DB)
  }
  
}

length(unique(ex_gwas1))
# 163
length(unique(ex_gwas2))
# 113
length(unique(ex_gwas3))
# 280
length(unique(ex_gwas4))
# 144
length(unique(ex_gwas5))
# 1
length(unique(ex_gem))
# 226

save(ex_gwas1, ex_gwas2, ex_gwas3, ex_gwas4, ex_gwas5, ex_gem, 
     file = "ex_cross_cohort_dups.RData")

rm(list = setdiff(ls(), c(ls(pattern="^cohort_"), ls(pattern="^ex_"))))

###########################################
### Combine exclude list of each cohort ###
###########################################

### gwas1
gwas1_dups <- data.frame(Sample = unique(ex_gwas1), 
                         Reason_for_Removal = "cross cohort duplicates")
gwas1_dups <- merge(gwas1_dups, cohort_dups[, c("AlzID_DB", "LabID_DB")], 
                    by.x = "Sample", by.y = "LabID_DB", all.x = TRUE)
gwas1_dups <- gwas1_dups[!duplicated(gwas1_dups), ]
gwas1_rm <- gwas1_dups
dim(gwas1_rm)
# [1] 163   3
table(gwas1_rm$Reason_for_Removal)
# cross cohort duplicates 
#                     163 

### gwas2 
gwas2_dups <- data.frame(Sample = unique(ex_gwas2), 
                         Reason_for_Removal = "cross cohort duplicates")
gwas2_dups <- merge(gwas2_dups, cohort_dups[, c("AlzID_DB", "LabID_DB")], 
                    by.x = "Sample", by.y = "LabID_DB", all.x = TRUE)
gwas2_dups <- gwas2_dups[!duplicated(gwas2_dups), ]

gwas2_low_quals <- openxlsx::read.xlsx("GWAS_Harmonization_exclu_Samples.xlsx", sheet = "GWAS2")
colnames(gwas2_dups) <- colnames(gwas2_low_quals)

gwas2_rm <- rbind(gwas2_dups, gwas2_low_quals)
gwas2_rm <- gwas2_rm[!duplicated(gwas2_rm), ]
dim(gwas2_rm)
# [1] 132   3
table(gwas2_rm$Reason_for_Removal)
#  cross cohort duplicates    identified duplicates           mislabeled sex missing call rate > 0.05 
#                      113                        7                       10                        2

### gwas3
gwas3_dups <- data.frame(Sample = unique(ex_gwas3), 
                         Reason_for_Removal = "cross cohort duplicates")
gwas3_dups <- merge(gwas3_dups, cohort_dups[, c("AlzID_DB", "LabID_DB")], 
                    by.x = "Sample", by.y = "LabID_DB", all.x = TRUE)
gwas3_dups <- gwas3_dups[!duplicated(gwas3_dups), ]

gwas3_low_quals <- openxlsx::read.xlsx("GWAS_Harmonization_exclu_Samples.xlsx", sheet = "GWAS3")
colnames(gwas3_dups) <- colnames(gwas3_low_quals)

gwas3_rm <- rbind(gwas3_dups, gwas3_low_quals)
gwas3_rm <- gwas3_rm[!duplicated(gwas3_rm), ]
dim(gwas3_rm)
# [1] 307   3
table(gwas3_rm$Reason_for_Removal)
#  cross cohort duplicates    identified duplicates missing call rate > 0.05 
#                      280                        2                       25 

### gwas4
gwas4_dups <- data.frame(Sample = unique(ex_gwas4), 
                         Reason_for_Removal = "cross cohort duplicates")
gwas4_dups <- merge(gwas4_dups, cohort_dups[, c("AlzID_DB", "LabID_DB")], 
                    by.x = "Sample", by.y = "LabID_DB", all.x = TRUE)
gwas4_dups <- gwas4_dups[!duplicated(gwas4_dups), ]

gwas4_low_quals <- openxlsx::read.xlsx("GWAS_Harmonization_exclu_Samples.xlsx", sheet = "GWAS4")
colnames(gwas4_dups) <- colnames(gwas4_low_quals)

gwas4_rm <- rbind(gwas4_dups, gwas4_low_quals)
gwas4_rm <- gwas4_rm[!duplicated(gwas4_rm), ]
dim(gwas4_rm)
# [1] 474   3
table(gwas4_rm$Reason_for_Removal)
# cross cohort duplicates    identified duplicates           mislabeled sex missing call rate > 0.05 
#                     144                       88                      106                      136

### gwas5
gwas5_dups <- data.frame(Sample = unique(ex_gwas5), 
                         Reason_for_Removal = "cross cohort duplicates")
gwas5_dups <- merge(gwas5_dups, cohort_dups[, c("AlzID_DB", "LabID_DB")], 
                    by.x = "Sample", by.y = "LabID_DB", all.x = TRUE)
gwas5_dups <- gwas5_dups[!duplicated(gwas5_dups), ]

gwas5_low_quals <- openxlsx::read.xlsx("GWAS_Harmonization_exclu_Samples.xlsx", sheet = "GWAS5")
colnames(gwas5_dups) <- colnames(gwas5_low_quals)

gwas5_rm <- rbind(gwas5_dups, gwas5_low_quals)
gwas5_rm <- gwas5_rm[!duplicated(gwas5_rm), ]
dim(gwas5_rm)
# [1] 15   3
table(gwas5_rm$Reason_for_Removal)
# cross cohort duplicates   identified duplicates          mislabeled sex 
#                       1                       3                      11

### gem
gem_dups <- data.frame(Sample = unique(ex_gem), 
                       Reason_for_Removal = "cross cohort duplicates")
gem_dups <- merge(gem_dups, cohort_dups[, c("AlzID_DB", "LabID_DB")], 
                  by.x = "Sample", by.y = "LabID_DB", all.x = TRUE)
gem_dups <- gem_dups[!duplicated(gem_dups), ]

gem_low_quals <- openxlsx::read.xlsx("GWAS_Harmonization_exclu_Samples.xlsx", sheet = "GEM")
colnames(gem_dups) <- colnames(gem_low_quals[1:3])

gem_rm <- rbind(gem_dups, gem_low_quals[,1:3])
gem_rm <- gem_rm[!duplicated(gem_rm), ]
dim(gem_rm)
# [1] 295   3
table(gem_rm$Reason_for_Removal)
# cross cohort duplicates    identified duplicates           mislabeled sex missing call rate > 0.05 
#                     226                       51                        7                       11 

library(xlsx)
write.xlsx(gwas1_rm, file="GWAS_Harmonization_exclu_Samples_Final_2024_01_08.xlsx", sheetName="GWAS1", row.names=FALSE)
write.xlsx(gwas2_rm, file="GWAS_Harmonization_exclu_Samples_Final_2024_01_08.xlsx", sheetName="GWAS2", append=TRUE, row.names=FALSE)
write.xlsx(gwas3_rm, file="GWAS_Harmonization_exclu_Samples_Final_2024_01_08.xlsx", sheetName="GWAS3", append=TRUE, row.names=FALSE)
write.xlsx(gwas4_rm, file="GWAS_Harmonization_exclu_Samples_Final_2024_01_08.xlsx", sheetName="GWAS4", append=TRUE, row.names=FALSE)
write.xlsx(gwas5_rm, file="GWAS_Harmonization_exclu_Samples_Final_2024_01_08.xlsx", sheetName="GWAS5", append=TRUE, row.names=FALSE)
write.xlsx(gem_rm, file="GWAS_Harmonization_exclu_Samples_Final_2024_01_08.xlsx", sheetName="GEM", append=TRUE, row.names=FALSE)

####################################################
### Modify exclude sample list based on VCF info ###
####################################################

gwas1_sample <- read.table("/zfs1/kfan/Ruyu/harmonization_Sep5/All_Cohorts/GWAS1_Sample.txt")
# fam <- read.table("/zfs1/kfan/Ruyu/harmonization_Sep5/GWAS_Clinical_Reference/GWAS_Reference_Data/Omni1-quad_07_01_2010_GenTrain_2_Kamboh_case_control.fam")
gwas1_ex_sample <- c()
for (i in 1:length(gwas1_rm$Sample)){
  print(i)
  gwas1_ex_sample <- c(gwas1_ex_sample, gwas1_sample[str_detect(gwas1_sample$V1, gwas1_rm$Sample[i]), ])
}
length(gwas1_rm$Sample) - length(gwas1_ex_sample)
# 55
write.table(as.data.frame(gwas1_ex_sample), file="GWAS1_ex_Samples.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

gwas2_sample <- read.table("/zfs1/kfan/Ruyu/harmonization_Sep5/All_Cohorts/GWAS2_Sample.txt")
gwas2_ex_sample <- c()
for (i in 1:length(gwas2_rm$Sample)){
  print(i)
  gwas2_ex_sample <- c(gwas2_ex_sample, gwas2_sample[str_detect(gwas2_sample$V1, gwas2_rm$Sample[i]), ])
}
length(gwas2_rm$Sample) - length(gwas2_ex_sample)
# 40
write.table(as.data.frame(gwas2_ex_sample), file="GWAS2_ex_Samples.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

gwas3_sample <- read.table("/zfs1/kfan/Ruyu/harmonization_Sep5/All_Cohorts/GWAS3_Sample.txt")
gwas3_ex_sample <- c()
for (i in 1:length(gwas3_rm$Sample)){
  print(i)
  gwas3_ex_sample <- c(gwas3_ex_sample, gwas3_sample[str_detect(gwas3_sample$V1, gwas3_rm$Sample[i]), ])
}
length(gwas3_rm$Sample) - length(gwas3_ex_sample)
# 209

gwas5_sample <- read.table("/zfs1/kfan/Ruyu/harmonization_Sep5/All_Cohorts/GWAS5_Sample.txt")
gwas5_ex_sample <- c()
for (i in 1:length(gwas5_rm$Sample)){
  print(i)
  gwas5_ex_sample <- c(gwas5_ex_sample, gwas5_sample[str_detect(gwas5_sample$V1, gwas5_rm$Sample[i]), ])
}
length(gwas5_rm$Sample) - length(gwas5_ex_sample)
# 4
write.table(as.data.frame(gwas5_ex_sample), file="GWAS5_ex_Samples.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

gem_sample <- read.table("/zfs1/kfan/Ruyu/harmonization_Sep5/All_Cohorts/GEM_Sample.txt")
gem_ex_sample <- c()
for (i in 1:length(gem_rm$Sample)){
  print(i)
  gem_ex_sample <- c(gem_ex_sample, gem_sample[str_detect(gem_sample$V1, gem_rm$Sample[i]), ])
}
length(gem_rm$Sample) - length(gem_ex_sample)
# 35
write.table(as.data.frame(gem_ex_sample), file="GEM_ex_Samples.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

save.image("cross_cohort_dups_check.RData")