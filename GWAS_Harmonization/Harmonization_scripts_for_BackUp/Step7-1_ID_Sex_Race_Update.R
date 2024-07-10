
### 2nd Imputation Sample Preparation

library(data.table)

setwd("/ix/kfan/Ruyu/harmonization/All_Cohorts")
################################################################################
### update FID/IID
################################################################################

batches = c(paste0("GWAS",1:5), "GEM")

fam <- fread("./AllCohorts_AllChr_Annot.fam")
dim(fam)

fam$V7 <- substr(fam$V1, 1, nchar(fam$V1)/2)

### read in batch lists

All_Cohorts_list <- list()
for (batch in batches){
  All_Cohorts_list[[batch]] <- read.table(paste0("/ix/kfan/Ruyu/harmonization/All_Cohorts/", batch, "_plink_IDs_noDup.txt"), sep = "\t", header = TRUE)
  All_Cohorts_list[[batch]]$Batch <- batch
}

All_Cohorts_dd <- do.call(rbind, All_Cohorts_list)
dim(All_Cohorts_dd)

table(All_Cohorts_dd$Batch)

# length(intersect(All_Cohorts_dd$V1, fam$V7))
All_Cohorts_2nd_Imputation <- merge(All_Cohorts_dd, fam[,c("V2","V7")], by.x = "V1", by.y = "V7")
dim(All_Cohorts_2nd_Imputation)

# during the conversion between VCF and PLINK format, the underscore in the ID column is causing some confusion, need to manually update the IDs for the final output

for (seq in 1:nrow(All_Cohorts_2nd_Imputation)){
  All_Cohorts_2nd_Imputation$V3[seq] <- substr(All_Cohorts_2nd_Imputation$V1[seq], 1, nchar(All_Cohorts_2nd_Imputation$V1[seq]) / 2)
}

for (seq in 1:nrow(All_Cohorts_2nd_Imputation)){
  
  scanName <- All_Cohorts_2nd_Imputation$V1[seq]
  
    All_Cohorts_2nd_Imputation$scanName[seq] <- substr(scanName, 1, nchar(scanName)/2)

}

write.table(All_Cohorts_2nd_Imputation[,c("scanName", "Batch", "V1", "V2")], "./All_Cohorts_2nd_Imputation_ID_Batch.txt", 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

ID_update_forPlink <- All_Cohorts_2nd_Imputation[,c("V2", "V2", "scanName", "scanName")]
row.names(ID_update_forPlink) <- NULL

write.table(ID_update_forPlink, 
            "./All_Cohorts_2nd_Imputation_ID_Update_forPlink.txt", 
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

################################################################################
### update Sex and Race
################################################################################

library(dplyr)

### read in an old DataRequest from Lynda
oldRef <- openxlsx::read.xlsx("/ix/kfan/Ruyu/harmonization/GWAS_Reference_Files/V_DataRequest_Ruyu_GWAS_Info_2023_12_20.xlsx")
dim(oldRef)
# [1] 11387    32

# head(oldRef)

check1 <- setdiff(All_Cohorts_2nd_Imputation$scanName, oldRef$Scan_Name)
# 41

# update oldRef ID
oldRef$Scan_Name2 <- gsub(" ", "", oldRef$Scan_Name)

check2 <- setdiff(All_Cohorts_2nd_Imputation$scanName, oldRef$Scan_Name2)
length(check2)
# 6

# update sex and race in oldRef
table(oldRef$Sex)
sum(is.na(oldRef$Sex))
table(oldRef$Sex_DB)
sum(is.na(oldRef$Sex_DB))

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
# Asian  Black Others  White 
#    27    710     74   9169
sum(is.na(oldRef$race_recode1))
# [1] 1407

# merge All_Cohorts_2nd_Imputation with oldRef 

d1 <- merge(All_Cohorts_2nd_Imputation[,c("Batch","V1","V2","scanName")], oldRef[,c("AlzID_DB", "LabID", "StudyID", "Scan_Name2", "sex_recode1", "race_recode1")], by.x = "scanName", by.y = "Scan_Name2") %>%
  distinct()
dim(d1)
# [1] 10944     9
colnames(d1) <- c("scanName","Batch","V1","V2","AlzID","LabID","StudyID","sex_recode","race_recode")

### read in the original GWAS clinical info
Info <- read.csv("/ix/kfan/Ruyu/harmonization/GWAS_Reference_Files/GWAS_ALL_2023_07_24.csv")
dim(Info)
# [1] 11113    16

check3 <- setdiff(check2, Info$LabID)

# recode sex and race
Info$sex_recode <- ifelse(Info$Sex < 0, NA, Info$Sex)
Info$race_recode <- ifelse(Info$Race %in% c("White", "W"), "White", 
                           ifelse(Info$Race %in% c("B", "Black"), "Black", 
                                  ifelse(Info$Race %in% c("A"), "Asian",
                                         ifelse(Info$Race < 0, NA, "Others"))))
# merge the missed ones in the dataRequest

d2 <- merge(All_Cohorts_2nd_Imputation[All_Cohorts_2nd_Imputation$scanName %in% check2, c("Batch","V1","V2","scanName")], 
            Info[,c("AlzID","LabID","STudyID", "sex_recode", "race_recode")], by.x = "scanName", by.y = "LabID") %>%
  distinct()
dim(d2)
# [1] 6 6
d2$LabID <- d2$scanName
colnames(d2) <- c("scanName","Batch","V1","V2","AlzID","StudyID","sex_recode","race_recode","LabID")

### rbind two parts of the d
d <- rbind(d1, 
            d2[,c("scanName","Batch","V1","V2","AlzID","LabID","StudyID","sex_recode","race_recode")]) %>%
  filter(!is.na(scanName)) %>%
  group_by(scanName) %>%
  filter(!(is.na(AlzID) & n() > 1)) %>%
  filter(!(is.na(LabID) & n() > 1))

table(d$Batch)
#  GEM GWAS1 GWAS2 GWAS3 GWAS4 GWAS5 
# 2731  2395   619   493  4280   235 

# sex check
table(d$sex_recode)
# sum(is.na(d$sex_recode)) # 1

# d[is.na(d$sex_recode),]$scanName
# "A0342_A0342"
d[d$scanName == "A0342_A0342",]$sex_recode <- Info[Info$LabID == "A0342", ]$Sex
table(d$sex_recode)
#    F    M 
# 6077 4676

# race check
table(d$race_recode)
sum(is.na(d$race_recode))
# 1341

noRace_scanName <- d[is.na(d$race_recode),]$scanName
Info$LabID2 <- paste0(Info$LabID,"_",Info$LabID)
for (seq in 1:length(noRace_scanName)){
  
  ID <- noRace_scanName[seq]
  
  if(nrow(Info[Info$LabID2 == ID, ]) != 0){
    d[d$scanName == ID,]$race_recode <- ifelse(Info[Info$LabID2 == ID, ]$race_recode > 0 & !is.na(Info[Info$LabID2 == ID, ]$race_recode), 
                                            Info[Info$LabID2 == ID, ]$race_recode, NA)
  } else {
    d[d$scanName == ID,]$race_recode <- NA
  }
}

table(d$race_recode)
# Asian  Black Others  White 
#    26    653     74   8660
sum(is.na(d$race_recode))
# 1340

d_noRace <- d[is.na(d$race_recode),]
table(d_noRace$Batch)
#  GWAS4 GWAS5 
#   1339     1 

### read in GWAS4 updated clinical info
GWAS4_Info <- read.csv("/ix/kfan/Ruyu/harmonization/GWAS_Reference_Files/GWAS_4_Clinical_Info_updated_2023_12_13.csv")
GWAS4_Info$V2 <- gsub(" ", "", GWAS4_Info$V1)
table(GWAS4_Info$Race)

GWAS4_Info$race_recode <- ifelse(GWAS4_Info$Race %in% c("White", "W"), "White", 
                                 ifelse(GWAS4_Info$Race %in% c("B", "Black"), "Black", 
                                        ifelse(GWAS4_Info$Race %in% c("A"), "Asian",
                                               ifelse(GWAS4_Info$Race < 0, NA, "Others"))))
table(GWAS4_Info$race_recode)
sum(is.na(GWAS4_Info$race_recode))
# 1400

check4 <- setdiff(d_noRace$scanName, GWAS4_Info$V2)

for (seq in 1:nrow(d[is.na(d$race_recode),])){
  
  labID <- d[is.na(d$race_recode),]$V2[seq]
  if(nrow(GWAS4_Info[GWAS4_Info$V2 == labID, ]) != 0){
    d[d$V2 == labID,]$race_recode <- ifelse(GWAS4_Info[GWAS4_Info$V2 == labID, ]$race_recode > 0 & !is.na(GWAS4_Info[GWAS4_Info$V2 == labID, ]$race_recode), 
                                            GWAS4_Info[GWAS4_Info$V2 == labID, ]$race_recode, NA)
  } else {
    d[d$V2 == labID,]$race_recode <- NA
  }
}

d_noRace <- d[is.na(d$race_recode),]
table(d_noRace$Batch)
# GWAS4 GWAS5 
#  1339     1

### output updated file
str(d)
dim(d)
# [1] 10753     9
table(d$Batch)
#  GEM GWAS1 GWAS2 GWAS3 GWAS4 GWAS5 
# 2731  2395   619   493  4280   235

# assign another group to those missing race
d$race_recode <- ifelse(is.na(d$race_recode), "Unknown", d$race_recode)
table(d$race_recode)

length(intersect(d$scanName, ID_update_forPlink$scanName))

write.table(d[,c("scanName", "AlzID", "LabID", "StudyID", "Batch", "sex_recode", "race_recode")], "./All_Cohorts_2nd_Imputation_ID_Batch_with_Sex_Race.txt", 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

################################################################################
### update Sex for Plink .fam file
################################################################################

# expects a file with FIDs and IIDs in the first two columns, 
# and sex information (1 or M = male, 2 or F = female, 0 = missing) in the (n+2)th column

d_sex_forPlink <- d[,c("scanName", "scanName", "sex_recode")]
table(d_sex_forPlink$sex_recode)
#    F    M 
# 6077 4676 
sum(is.na(d_sex_forPlink$sex_recode))
# 0 

write.table(d_sex_forPlink, "./All_Cohorts_2nd_Imputation_Sex_Update_forPlink.txt", 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
