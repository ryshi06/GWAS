
### remove cross-cohort duplicated samples

batches = c(paste0("GWAS",1:5), "GEM")

scan_name <- list()
for (batch in batches){
  scan_name[[batch]] <- read.table(paste0("/ix/kfan/Ruyu/harmonization/All_Cohorts/pseudoGWAS_input/", batch,"_allChr_scanName.txt"))
  scan_name[[batch]]$Batch <- batch
}

scan_name_d <- do.call(rbind, scan_name)
table(scan_name_d$Batch)
#  GEM GWAS1 GWAS2 GWAS3 GWAS4 GWAS5 
# 2738  2395   711   613  4318   245 

length(unique(scan_name_d$V1))
# 10816

scan_name_dups <- scan_name_d[scan_name_d$V1 %in% scan_name_d$V1[duplicated(scan_name_d$V1)],]
# [1] 392   2
table(scan_name_dups$Batch)
# GEM GWAS2 GWAS3 GWAS4 
# 165    84   121    22 

scan_name_dups_keep <- data.frame(V1 = unique(scan_name_dups$V1),
                                  Batch = NA)
for (sample in scan_name_dups_keep$V1) {
  batch_list <- scan_name_dups[scan_name_dups$V1 == sample, "Batch"]
  if ("GWAS5" %in% batch_list) {
    scan_name_dups_keep[scan_name_dups_keep$V1 == sample, "Batch"] <- "GWAS5"
  } else if ("GEM" %in% batch_list) {
    scan_name_dups_keep[scan_name_dups_keep$V1 == sample, "Batch"] <- "GEM"
  } else if ("GWAS4" %in% batch_list) {
    scan_name_dups_keep[scan_name_dups_keep$V1 == sample, "Batch"] <- "GWAS4"
  } else if ("GWAS3" %in% batch_list) {
    scan_name_dups_keep[scan_name_dups_keep$V1 == sample, "Batch"] <- "GWAS3"
  } else if ("GWAS2" %in% batch_list) {
    scan_name_dups_keep[scan_name_dups_keep$V1 == sample, "Batch"] <- "GWAS2"
  } else {
    scan_name_dups_keep[scan_name_dups_keep$V1 == sample, "Batch"] <- "GWAS1"
  }
}

scan_name_noDup <- rbind(scan_name_d[!(scan_name_d$V1 %in% scan_name_d$V1[duplicated(scan_name_d$V1)]),],
                         scan_name_dups_keep)
table(scan_name_noDup$Batch)
#  GEM GWAS1 GWAS2 GWAS3 GWAS4 GWAS5 
# 2738  2395   627   493  4318   245 

for (batch in batches){
  write.table(scan_name_noDup[scan_name_noDup$Batch == batch,], 
              paste0("/ix/kfan/Ruyu/harmonization/All_Cohorts/", batch, "_plink_IDs_noDup.txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

### extract IDs of the duplicates
##GWAS2_dups <- setdiff(scan_name_d[scan_name_d$Batch == "GWAS2",]$V1, scan_name_noDup[scan_name_noDup$Batch == "GWAS2",]$V1)
##GWAS3_dups <- setdiff(scan_name_d[scan_name_d$Batch == "GWAS3",]$V1, scan_name_noDup[scan_name_noDup$Batch == "GWAS3",]$V1)

### remove sex-mismatch samples

setwd("/ix/kfan/Ruyu/harmonization/All_Cohorts")

batches = c(paste0("GWAS",1:5), "GEM")

excl.list <- list()

cross_batch_dups <- list()
sex_mismatch <- list()
high_missing_call_rate <- list()
identified_dups <- list()

## terms <- c("cross cohort duplicates", "identified duplicates", "mislabeled sex", "missing call rate > 0.05")

for (batch in batches){
  excl.list[[batch]] <- openxlsx::read.xlsx("GWAS_Harmonization_exclu_Samples_Final_2024_05_01.xlsx", sheet = batch)
  
  excl.list[[batch]]$batch <- batch
  colnames(excl.list[[batch]]) <- c("Sample", "Reason_for_Removal", "AlzID", "Batch")
  cross_batch_dups[[batch]] <- excl.list[[batch]][excl.list[[batch]]$Reason_for_Removal == "cross cohort duplicates",]
  sex_mismatch[[batch]] <- excl.list[[batch]][excl.list[[batch]]$Reason_for_Removal == "mislabeled sex",]
  identified_dups[[batch]] <- excl.list[[batch]][excl.list[[batch]]$Reason_for_Removal == "identified duplicates",]
  high_missing_call_rate[[batch]] <- excl.list[[batch]][excl.list[[batch]]$Reason_for_Removal == "missing call rate > 0.05",]
}

sex_mismatch_dd <- do.call(rbind, sex_mismatch)
sex_mismatch_dd$V1 <- paste0(sex_mismatch_dd$Sample, "_", sex_mismatch_dd$Sample)
sex_mismatch_list <- unique(sex_mismatch_dd$Sample)
sex_mismatch_rm <- data.frame(FID = paste0(sex_mismatch_list, "_", sex_mismatch_list),
                                 IID = paste0(sex_mismatch_list, "_", sex_mismatch_list))
write.table(sex_mismatch_rm, "All_Cohorts_ScanName_Sex_Mismatch_plink_ex.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

noDups.fam <- fread("./pseudoGWAS_input/All_Cohorts_allChr_noDups.fam")
setdiff(sex_mismatch_rm$FID, noDups.fam$V1)

# View(sex_mismatch_dd[sex_mismatch_dd$V1 %in% setdiff(sex_mismatch_rm$FID, noDups.fam$V1), ])
# View(sex_mismatch_dd[sex_mismatch_dd$V1 %in% intersect(sex_mismatch_rm$FID, noDups.fam$V1), ])

write.xlsx(sex_mismatch_dd[sex_mismatch_dd$V1 %in% intersect(sex_mismatch_rm$FID, noDups.fam$V1), ], 
           "./Harmonization_Sample_Ex_Sex_Mismatch.xlsx", colNames = TRUE)


##table(scan_name_d$Batch)
# GEM GWAS1 GWAS2 GWAS3 GWAS4 GWAS5 
# 2738  2395   711   613  4318   245 
##table(scan_name_noDup$Batch)
#  GEM GWAS1 GWAS2 GWAS3 GWAS4 GWAS5 
# 2738  2395   627   493  4318   245 


##GWAS2_RM <- scan_name_d[scan_name_d$V1 %in% setdiff(scan_name_d[scan_name_d$Batch == "GWAS2",]$V1, 
##                    scan_name_noDup[scan_name_noDup$Batch == "GWAS2",]$V1),]

##table(GWAS2_RM$Batch)
# GEM GWAS2 GWAS3 
# 83    84    17

##GWAS3_RM <- scan_name_d[scan_name_d$V1 %in% setdiff(scan_name_d[scan_name_d$Batch == "GWAS3",]$V1, 
##                                                    scan_name_noDup[scan_name_noDup$Batch == "GWAS3",]$V1),]

##table(GWAS3_RM$Batch)
# GEM GWAS2 GWAS3 GWAS4 
# 98    16   120    22

##GWAS3_RM_gwas2 <- GWAS3_RM[GWAS3_RM$Batch == "GWAS2",]
##table(GWAS2_RM$Batch)
##table(GWAS2_RM[!(GWAS2_RM$V1 %in% GWAS3_RM_gwas2$V1), ]$Batch)

##table(GWAS3_RM$Batch)
##table(GWAS3_RM[!(GWAS3_RM$V1 %in% GWAS3_RM_gwas2$V1), ]$Batch)

# Note: those 16 were genotyped three times (GEM, GWAS2, and GWAS3)