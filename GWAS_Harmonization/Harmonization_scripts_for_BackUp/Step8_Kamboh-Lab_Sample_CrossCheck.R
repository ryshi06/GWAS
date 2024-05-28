
rm(list = ls())
setwd("/Users/shiruyu/Desktop/tmp")

library(data.table)
library(dplyr)

### Note: this script is used to check the data consistencies among GWAS Harmonization Table, APOE24 Table, and Autopsy 

###--------------------------------------------------------------------------###
### load required files
###--------------------------------------------------------------------------###

### the original file sent from Lynda related to harmonization
oldRef <- openxlsx::read.xlsx("V_DataRequest_Ruyu_GWAS_Info_2023_12_20.xlsx")
length(unique(oldRef$Scan_Name))

oldRef$scanName <- gsub(" ", "", oldRef$Scan_Name)
oldRef$sex_recode1 <- ifelse(oldRef$Sex_DB < 0, NA, oldRef$Sex_DB)
oldRef$sex_recode <- ifelse(!is.na(oldRef$sex_recode1), oldRef$sex_recode1, oldRef$Sex)
oldRef$race_recode1 <- ifelse(!is.na(oldRef$Race_DB), oldRef$Race_DB, oldRef$Race)
oldRef$race_recode <- ifelse(oldRef$race_recode1 < 0 | is.na(oldRef$race_recode1), "Unknown",
                             ifelse(oldRef$race_recode1 %in% c("W", "White"), "White", 
                                    ifelse(oldRef$race_recode1 %in% c("B", "Black"), "Black",
                                           ifelse(oldRef$race_recode1 %in% c("A", "Asian"), "Asian","Others"))))
oldRef$apoe_recode1 <- ifelse(oldRef$APOE_DB < 0 | oldRef$APOE_DB == 99 | is.na(oldRef$APOE_DB), NA, oldRef$APOE_DB)
oldRef$apoe_recode <- ifelse(!is.na(oldRef$apoe_recode1), oldRef$apoe_recode1, oldRef$APOE)
oldRef$case_control1 <- ifelse(oldRef$Case_Control_DB == "-1" | is.na(oldRef$APOE_DB), NA, oldRef$Case_Control_DB)
oldRef$case_control <-  ifelse(!is.na(oldRef$case_control1), oldRef$case_control1, oldRef$Case_Control)

### the APOE24 table
apoe24 <- openxlsx::read.xlsx("APOE24_Sample_Characteristics_for_TabPrep_N15186.xlsx")
apoe24_split <- apoe24 %>%
  tidyr::separate_rows(LabIDList, sep = ",") %>%
  as.data.frame()

apoe24_split$labID_numbers <- gsub("[^0-9]", "", apoe24_split$LabIDList)

apoe24_split$race_recode <- ifelse(apoe24_split$Race < 0 | is.na(apoe24_split$Race), "Unknown",
                             ifelse(apoe24_split$Race %in% c("W", "White"), "White", 
                                    ifelse(apoe24_split$Race %in% c("B", "Black"), "Black",
                                           ifelse(apoe24_split$Race %in% c("A", "Asian"), "Asian","Others"))))
apoe24_split$sex_recode <- ifelse(apoe24_split$Sex < 0, NA, apoe24_split$Sex)
apoe24_split$case_control <- ifelse(apoe24_split$Case_Control_Status == "-1", NA, apoe24_split$Case_Control_Status)

### the 2nd imputation phenotype file
pheno <- fread("All_Cohorts_2nd_Imputation_Phenotype.txt")

### the final harmonized phenotype file
pheno.final  <- fread("AllCohorts_AllChr_Annot_QC5.fam")
# View(pheno.final)

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

gwas1.fam <- fread("/Users/shiruyu/Desktop/tmp/Omni1-quad_07_01_2010_GenTrain_2_Kamboh_case_control.fam")
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

### read in the ibd removal info
ibd_rm <- fread("IBD_Remove_list.txt")

### get the total N = 11387 list

all <- rbind(ex2[,c("V1", "Batch")], allScan)
table(all$Batch)
#  GEM GWAS1 GWAS2 GWAS3 GWAS4 GWAS5 
# 2800  2440   720   640  4539   248
# View(all)

priority_order <- c("GWAS5", "GEM", "GWAS4", "GWAS3", "GWAS2", "GWAS1")

###--------------------------------------------------------------------------###
### create a new dataframe to summarize all related info
###--------------------------------------------------------------------------###

pheno.final$QC_Check <- "Pass"
dim(pheno.final)

dd <- merge(pheno, pheno.final[,c("V1", "QC_Check")], by.x = "scanName", by.y = "V1", all.x = TRUE)
dd[is.na(dd$QC_Check),]$QC_Check <- "Failed"
table(dd$QC_Check) 
# Failed   Pass 
#    142  10611

dd$`Removal.Reasons` <- ifelse(dd$QC_Check == "Failed", "Cross-batch Duplicates by IBD > 0.9", NA)

dim(sex.ex)
# add back sex mismatched samples
sex.dd <- merge(sex.ex[,c("SampleID","Removal.Reasons","Batch")], oldRef[,c("AlzID","LabID","StudyID","sex_recode","race_recode","apoe_recode", "scanName","case_control")], by.x = "SampleID", by.y = "scanName", all.x = TRUE) %>% distinct()
sex.dd$SampleID[duplicated(sex.dd$SampleID)]
sex.dd_filtered <- sex.dd %>%
  group_by(SampleID) %>%
  # Check if V1 is duplicated
  mutate(duplicated_V1 = n() > 1) %>%
  # Remove rows with missing values in AlzID, LabID, or StudyID for duplicated V1s
  filter(!(duplicated_V1 & (is.na(AlzID) | is.na(LabID) | is.na(StudyID)))) %>%
  # Select the columns you want to keep
  select(-duplicated_V1) %>%
  as.data.frame()

sex.dd_filtered[sex.dd_filtered$SampleID == "4210301",]$sex_recode <- apoe24_split[apoe24_split$LabIDList == "4210301",]$sex_recode
sex.dd_filtered[sex.dd_filtered$SampleID == "4210301",]$race_recode <- apoe24_split[apoe24_split$LabIDList == "4210301",]$race_recode
sex.dd_filtered[sex.dd_filtered$SampleID == "4210301",]$apoe_recode <- apoe24_split[apoe24_split$LabIDList == "4210301",]$APOE
sex.dd_filtered[sex.dd_filtered$SampleID == "4210301",]$case_control <- apoe24_split[apoe24_split$LabIDList == "4210301",]$Case_Control_Status

sex.dd_filtered[sex.dd_filtered$SampleID == "5206154",]$sex_recode <- apoe24_split[apoe24_split$LabIDList == "5206154",]$sex_recode
sex.dd_filtered[sex.dd_filtered$SampleID == "5206154",]$race_recode <- apoe24_split[apoe24_split$LabIDList == "5206154",]$race_recode
sex.dd_filtered[sex.dd_filtered$SampleID == "5206154",]$apoe_recode <- apoe24_split[apoe24_split$LabIDList == "5206154",]$APOE
sex.dd_filtered[sex.dd_filtered$SampleID == "5206154",]$case_control <- apoe24_split[apoe24_split$LabIDList == "5206154",]$Case_Control_Status
sex.dd_filtered$QC_Check <- "Failed"

colnames(sex.dd_filtered) <- c("scanName", "Removal.Reasons", "Batch","AlzID", "LabID","StudyID", "sex_recode", "race_recode", "apoe_recode", "case_control","QC_Check")
dd1 <- rbind(dd[,c("scanName","AlzID", "LabID", "StudyID", "Batch", "QC_Check", "Removal.Reasons")], 
             sex.dd_filtered[,c("scanName","AlzID", "LabID", "StudyID", "Batch", "QC_Check", "Removal.Reasons")])
table(dd1$QC_Check)

allScan.dd <- allScan %>%
  # Create a priority rank based on the priority_order
  mutate(priority_rank = match(Batch, priority_order)) %>%
  # Group by ID to process each group separately
  group_by(V1) %>%
  # Arrange by priority_rank within each group
  arrange(priority_rank) %>%
  # Mark the row with the highest priority (lowest priority_rank) within each group
  mutate(QC_Check = if_else(row_number() == 1, "Pass", "Failed")) %>%
  # Ungroup to return to regular data frame operations
  ungroup() %>%
  # Select the final columns needed
  select(V1, Batch, QC_Check)

allScan.dd$"Removal.Reasons" <- ifelse(allScan.dd$QC_Check == "Failed", "Cross-batch Duplicates by overlapped Sample ID across batches", NA)
allScan.dd1 <- merge(allScan.dd, oldRef[!is.na(oldRef$AlzID),c("AlzID","LabID","StudyID","scanName")], by.x = "V1", by.y = "scanName", all.x = TRUE) %>%
  distinct()

allScan.dd1_filtered <- allScan.dd1 %>%
  group_by(V1) %>%
  # Check if V1 is duplicated
  mutate(duplicated_V1 = n() > 1) %>%
  # Remove rows with missing values in AlzID, LabID, or StudyID for duplicated V1s
  filter(!(duplicated_V1 & (is.na(AlzID) | is.na(LabID) | is.na(StudyID)))) %>%
  # Select the columns you want to keep
  select(-duplicated_V1) %>%
  as.data.frame()

head(allScan.dd1_filtered)
colnames(allScan.dd1_filtered) <- c("scanName", "Batch","QC_Check","Removal.Reasons","AlzID","LabID","StudyID")

dd2 <- rbind(dd1, allScan.dd1_filtered[allScan.dd1_filtered$QC_Check == "Failed", colnames(allScan.dd1_filtered) %in% colnames(dd1)])
table(dd2$QC_Check)

ex2$QC_Check <- "Failed"

# length(intersect(ex2$V1, oldRef$scanName))
ex2$V2 <- ifelse(ex2$Batch == "GWAS1", sub("^[^_]*_", "", ex2$V1), ex2$V1)
# setdiff(ex2$V2, oldRef$scanName)

ex3 <- merge(ex2, oldRef[,c("AlzID","LabID","StudyID","scanName")], by.x = "V2", by.y = "scanName") %>% distinct()
ex3$V1[duplicated(ex3$V1)]

ex3[ex3$V1 == "A0035_A0035",]$LabID <- "A0035"
ex3[ex3$V1 == "A0035",]$LabID <- "A0035"
ex3[ex3$V1 == "A0630_A0630",]$LabID <- "A0630"
ex3[ex3$V1 == "A0630",]$LabID <- "A0630"
ex3[ex3$V1 == "A0225_A0225",]$LabID <- "A0225"
ex3[ex3$V1 == "A0824_A0824",]$LabID <- "A0824"
ex3[ex3$V1 == "A2167_A2167",]$LabID <- "A2167"

ex3_filtered <- ex3 %>%
  group_by(V1) %>%
  # Check if V1 is duplicated
  mutate(duplicated_V1 = n() > 1) %>%
  # Remove rows with missing values in AlzID, LabID, or StudyID for duplicated V1s
  filter(!(duplicated_V1 & (is.na(AlzID) | is.na(LabID) | is.na(StudyID)))) %>%
  # Select the columns you want to keep
  select(-duplicated_V1) %>%
  as.data.frame()

# head(ex3_filtered)
colnames(ex3_filtered) <- c("V2", "scanName", "Removal.Reasons", "Batch", "QC_Check","AlzID", "LabID","StudyID")
# head(dd2)

dd3 <- rbind(dd2, ex3_filtered[,colnames(dd2)])
table(dd3$Batch)
#  GEM GWAS1 GWAS2 GWAS3 GWAS4 GWAS5 
# 2800  2440   720   640  4539   248
table(dd3$Removal.Reasons)

###--------------------------------------------------------------------------###
### combine the info with the references (add phenotype)
###--------------------------------------------------------------------------###

d_11387 <- dd3
colnames(d_11387) <- c("scanName","AlzID","LabID","StudyID","Batch","GWAS_QC","GWAS_QC_Failed_Notes")
head(d_11387)

# update race in oldRef
table(oldRef$Study_DB)
oldRef$race_recode <- ifelse(oldRef$Study_DB == "BENIN", "Black", oldRef$race_recode)
oldRef$race_recode <- ifelse(oldRef$Study_DB == "SLVDS", "White", oldRef$race_recode)

# update race in apoe24 table
apoe24_split$race_recode <- ifelse(apoe24_split$StudyList == "BENIN", "Black", apoe24_split$race_recode)
apoe24_split$race_recode <- ifelse(apoe24_split$StudyList == "SLVDS", "White", apoe24_split$race_recode)

d_11387$Study <- NA
d_11387$Sex <- NA
d_11387$Race <- NA
d_11387$APOE <- NA
d_11387$Case_Control_Status <- NA

d_11387 <- d_11387 %>%
  dplyr::mutate(across(c(Study, Sex, Race, APOE, Case_Control_Status), as.character))

for (seq in 1:nrow(d_11387)){
  
  ID <- d_11387[seq,]$scanName
  print(paste0("Processing row ", seq, " out of ", nrow(d_11387),"."))
  
  tmp1 <- pheno[pheno$scanName == ID,]
  
  if(nrow(tmp1) > 0){
    
    d_11387[d_11387$scanName == ID,]$LabID <- ifelse(is.na(d_11387[d_11387$scanName == ID,]$LabID), tmp1$LabID, d_11387[d_11387$scanName == ID,]$LabID)
    d_11387[d_11387$scanName == ID,]$StudyID <- ifelse(is.na(d_11387[d_11387$scanName == ID,]$StudyID), tmp1$StudyID, d_11387[d_11387$scanName == ID,]$StudyID)
    
    d_11387[d_11387$scanName == ID,]$Sex <- tmp1$sex_recode
    d_11387[d_11387$scanName == ID,]$Race <- ifelse(tmp1$race_recode != "Unknown", tmp1$race_recode, NA)
    d_11387[d_11387$scanName == ID,]$APOE <- tmp1$apoe_recode
    d_11387[d_11387$scanName == ID,]$Case_Control_Status <- tmp1$case_control
  }
  
  tmp2 <- apoe24_split[apoe24_split$LabIDList %in% c(toupper(ID), tolower(ID)),]
  if(nrow(tmp2) > 0){
    
    d_11387[d_11387$scanName == ID,]$LabID <- ifelse(is.na(d_11387[d_11387$scanName == ID,]$LabID), tmp2$LabIDList, d_11387[d_11387$scanName == ID,]$LabID)
    
    d_11387[d_11387$scanName == ID,]$Study <- ifelse(is.na(d_11387[d_11387$scanName == ID,]$Study), tmp2$StudyList, d_11387[d_11387$scanName == ID,]$Study)
    d_11387[d_11387$scanName == ID,]$Sex <- ifelse(is.na(d_11387[d_11387$scanName == ID,]$Sex), tmp2$sex_recode, d_11387[d_11387$scanName == ID,]$Sex)
    d_11387[d_11387$scanName == ID,]$Race <- ifelse(is.na(d_11387[d_11387$scanName == ID,]$Race), tmp2$race_recode, d_11387[d_11387$scanName == ID,]$Race)
    d_11387[d_11387$scanName == ID,]$APOE <- ifelse(is.na(d_11387[d_11387$scanName == ID,]$APOE), tmp2$APOE, d_11387[d_11387$scanName == ID,]$APOE)
    d_11387[d_11387$scanName == ID,]$Case_Control_Status <- ifelse(is.na(d_11387[d_11387$scanName == ID,]$Case_Control_Status), tmp2$case_control, d_11387[d_11387$scanName == ID,]$Case_Control_Status)
  }
  
  tmp3 <- oldRef[oldRef$scanName == ID, c("scanName","AlzID","LabID","StudyID","Study_DB","sex_recode","race_recode")] %>% distinct()
  if(nrow(tmp3) > 0){
    
    d_11387[d_11387$scanName == ID,]$LabID <- ifelse(is.na(d_11387[d_11387$scanName == ID,]$LabID), unique(tmp3$LabID[!is.na(tmp3$LabID)]), d_11387[d_11387$scanName == ID,]$LabID)
    d_11387[d_11387$scanName == ID,]$StudyID <- ifelse(is.na(d_11387[d_11387$scanName == ID,]$StudyID), unique(tmp3$StudyID[!is.na(tmp3$StudyID)]), d_11387[d_11387$scanName == ID,]$StudyID)
    d_11387[d_11387$scanName == ID,]$Study <- ifelse(length(unique(tmp3$Study_DB[!is.na(tmp3$Study_DB)]))>0, unique(tmp3$Study_DB[!is.na(tmp3$Study_DB)]), NA)
    d_11387[d_11387$scanName == ID,]$Sex <- ifelse(is.na(d_11387[d_11387$scanName == ID,]$Sex), tmp3$sex_recode, d_11387[d_11387$scanName == ID,]$Sex)
    d_11387[d_11387$scanName == ID,]$Race <- ifelse(is.na(d_11387[d_11387$scanName == ID,]$Race), tmp3$race_recode, d_11387[d_11387$scanName == ID,]$Race)
    
  }
}

### Study
sum(is.na(d_11387$Study))
# d_11387[is.na(d_11387$Study),]

d_11387[d_11387$LabID %in% c("M0279","M0280"),]$Study <- "MOVIES"
d_11387[d_11387$LabID == "6208185",]$Study <- "GEM"
d_11387[d_11387$LabID == "G0444",]$Study <- "MYHAT"
d_11387[d_11387$LabID %in% d_11387[is.na(d_11387$Study),]$LabID[grep("A|B",d_11387[is.na(d_11387$Study),]$LabID)],]$Study <- "ADRC"

d_11387[is.na(d_11387$Study),]$scanName
# [1] "S275" "S313"

table(d_11387$Study)
# ADRC      BACNE      BENIN        GEM        HCP Heartscore     IGNITE       LLDP     MOVIES      MYHAT       SARP 
# 3296        158        787       3016        230        150        673         26        777       1574          7 
# SENIORS      SLVDS      WHITE 
#      21        576         94 

### Sex
d_11387$Sex <- ifelse(d_11387$Sex < 0, NA, d_11387$Sex)
sum(is.na(d_11387$Sex))

dd_11387 <- d_11387 %>%
  group_by(LabID) %>%
  mutate(
    Sex = if_else(is.na(LabID), Sex, coalesce(Sex, first(Sex[!is.na(Sex)]))),
    Race = if_else(is.na(LabID), Race, coalesce(Race, first(Race[!is.na(Race)]))),
    APOE = if_else(is.na(LabID), APOE, coalesce(APOE, first(APOE[!is.na(APOE)]))),
    Case_Control_Status = if_else(is.na(LabID), Case_Control_Status, coalesce(Case_Control_Status, first(Case_Control_Status[!is.na(Case_Control_Status)])))
  ) %>%
  ungroup()

table(dd_11387$Sex)
#    F    M 
# 6387 4998

### Race
dd_11387$Race <- ifelse(dd_11387$Race < 0 | dd_11387$Race == "Unknown", NA, dd_11387$Race)
sum(is.na(dd_11387$Race)) # 41

table(dd_11387[is.na(dd_11387$Race),]$Batch)
# GWAS1 GWAS4 GWAS5 
#     2    35     4
table(dd_11387[is.na(dd_11387$Race),]$Study)
# BACNE        HCP Heartscore     IGNITE     MOVIES      MYHAT    SENIORS 
#     1          1          1          9          2          4         21

sum(is.na(dd_11387[is.na(dd_11387$Race),]$Study))
# 2

table(dd_11387$Race)
# Asian  Black Others  White 
#    27   1498     74   9747

### APOE
unique(dd_11387$APOE)

table(dd_11387$APOE)
# 22   23   24   33   34   44 
# 53 1118  268 6006 3015  483

# sum(is.na(dd_11387$APOE)) 
# 444
# dd_11387[is.na(dd_11387$APOE),]

### Case-Control 
table(dd_11387$Case_Control_Status)
dd_11387$Case_Control_Status[dd_11387$Case_Control_Status == "AD"] <- "Case"
table(dd_11387$Case_Control_Status)
# Case         Control Dementia-Non-AD 
# 3405            7847               7 
sum(is.na(dd_11387$Case_Control_Status))
table(dd_11387[is.na(dd_11387$Case_Control_Status),]$Batch)
# GWAS1 GWAS2 GWAS3 GWAS4 GWAS5 
#    23     2     4    91     8  
dd_11387[is.na(dd_11387$Case_Control_Status),]$scanName

dd_11387$GWAS <- ifelse(dd_11387$GWAS_QC == "Pass", "Y", "N")
table(dd_11387$GWAS)
#   N     Y 
# 776 10611

openxlsx::write.xlsx(dd_11387[,c("AlzID", "scanName", "LabID", "StudyID", "Study", "Batch", "GWAS", "GWAS_QC", "GWAS_QC_Failed_Notes", "Sex", "Race", "APOE", "Case_Control_Status")], 
                     "./GWAS_Harmonization_N=11387_Phenotype_2024_05_26.xlsx")

###--------------------------------------------------------------------------###
### combine the info with the APOE24 table
###--------------------------------------------------------------------------###

# head(apoe24_split)
# head(dd_11387)

apoe24_split$Data_Source <- "APOE24_Table_N15816"

length(intersect(dd_11387$AlzID, apoe24_split$AlzID))
# 10470

# merge by AlzID
d1 <- merge(dd_11387[dd_11387$AlzID %in% intersect(dd_11387$AlzID, apoe24_split$AlzID),], 
            apoe24_split[apoe24_split$AlzID %in% intersect(dd_11387$AlzID, apoe24_split$AlzID), c("AlzID", "StudyList","LabIDList","Data_Source")], 
            by = "AlzID") %>% distinct()
length(unique(d1$AlzID))
# 10470 
length(unique(d1$scanName))
# 10708
d1$AlzID.Harmonization <- d1$AlzID
d1$AlzID.APOE24 <- d1$AlzID
colnames(d1) <- c("AlzID","scanName","LabID.Harmonization","StudyID","Batch","GWAS_QC","GWAS_QC_Failed_Notes","Study.Harmonization","Sex","Race",
                  "APOE","Case_Control_Status","GWAS","Study.APOE24","LabID.APOE24","Data_Source","AlzID.Harmonization","AlzID.APOE24")

# for those with unmatched AlzID, check LabID
diff1 <- setdiff(dd_11387$AlzID, apoe24_split$AlzID)
length(diff1) # 260

diff1_LabIDs <- dd_11387[dd_11387$AlzID %in% diff1, ]$LabID[!is.na(dd_11387[dd_11387$AlzID %in% diff1, ]$LabID)]
# 193
d2 <- merge(dd_11387[dd_11387$LabID %in% diff1_LabIDs, ],
            apoe24_split[apoe24_split$LabIDList %in% diff1_LabIDs, c("AlzID", "StudyList","LabIDList","Data_Source")],
            by.x = "LabID", by.y = "LabIDList") %>% distinct()
colnames(d2) <- c("LabID","scanName","AlzID.Harmonization","StudyID","Batch","GWAS_QC","GWAS_QC_Failed_Notes","Study.Harmonization","Sex","Race","APOE","Case_Control_Status","GWAS","AlzID.APOE24", "Study.APOE24", "Data_Source")
length(unique(d2$scanName))
# 212
d2$LabID.Harmonization <- d2$LabID
d2$LabID.APOE24 <- d2$LabID

# for those scanNames not in d1 and d2, check StudyID
scanName_diff <- unique(setdiff(unique(dd_11387$scanName), unique(c(d1$scanName, d2$scanName))))
# 281
# length(intersect(toupper(scanName_diff), toupper(apoe24_split$LabIDList))) # 0
# length(intersect(toupper(dd_11387[dd_11387$scanName %in% scanName_diff, ]$StudyID), toupper(apoe24_split$LabIDList))) # 0
d3 <- dd_11387[dd_11387$scanName %in% scanName_diff, ]
colnames(d3) <- c("scanName", "AlzID.Harmonization", "LabID.Harmonization", "StudyID", "Batch", "GWAS_QC", "GWAS_QC_Failed_Notes","Study.Harmonization", "Sex","Race","APOE","Case_Control_Status","GWAS")
d3$AlzID.APOE24 <- NA
d3$Data_Source <- NA
d3$Study.APOE24 <- NA
d3$LabID.APOE24 <- NA

# for those LabIDLists not in d1, d2, and d3 and not overlapped with scanName
diff2_LabIDs <- setdiff(unique(apoe24_split$LabIDList), unique(c(d1$LabID.APOE24, d2$LabID.APOE24, d3$LabID.APOE24, dd_11387$scanName)))
d4 <- apoe24_split[apoe24_split$LabIDList %in% diff2_LabIDs, c("AlzID", "StudyList","LabIDList","Data_Source")]
colnames(d4) <- c("AlzID.APOE24", "Study.APOE24", "LabID.APOE24", "Data_Source")

d4$scanName <- NA
d4$AlzID.Harmonization <- NA
d4$LabID.Harmonization <- NA
d4$StudyID <- NA
d4$Batch <- NA
d4$GWAS <- "Not Available"
d4$GWAS_QC <- NA
d4$GWAS_QC_Failed_Notes <- NA
d4$Study.Harmonization <- NA
d4$Sex <- NA
d4$Race <- NA
d4$APOE <- NA
d4$Case_Control_Status <- NA

d <- rbind(d1[,c("scanName","AlzID.Harmonization","LabID.Harmonization","StudyID","Batch","Study.Harmonization",
                 "GWAS","GWAS_QC","GWAS_QC_Failed_Notes","Sex","Race","APOE","Case_Control_Status",
                 "AlzID.APOE24","LabID.APOE24","Study.APOE24","Data_Source")],
           d2[,c("scanName","AlzID.Harmonization","LabID.Harmonization","StudyID","Batch","Study.Harmonization",
                 "GWAS","GWAS_QC","GWAS_QC_Failed_Notes","Sex","Race","APOE","Case_Control_Status",
                 "AlzID.APOE24","LabID.APOE24","Study.APOE24","Data_Source")],
           d3[,c("scanName","AlzID.Harmonization","LabID.Harmonization","StudyID","Batch","Study.Harmonization",
                 "GWAS","GWAS_QC","GWAS_QC_Failed_Notes","Sex","Race","APOE","Case_Control_Status",
                 "AlzID.APOE24","LabID.APOE24","Study.APOE24","Data_Source")],
           d4[,c("scanName","AlzID.Harmonization","LabID.Harmonization","StudyID","Batch","Study.Harmonization",
                 "GWAS","GWAS_QC","GWAS_QC_Failed_Notes","Sex","Race","APOE","Case_Control_Status",
                 "AlzID.APOE24","LabID.APOE24","Study.APOE24","Data_Source")]) %>% distinct() 
dim(d)
# [1] 17108    17
d$N15816_APOE24 <- ifelse(!is.na(d$Data_Source), "Y", "N")
table(d$N15816_APOE24)
#   N     Y 
# 282 16826

length(unique(intersect(d$AlzID.APOE24, apoe24_split$AlzID)))
# 15064
length(setdiff(apoe24_split$AlzID, d$AlzID.APOE24))
# 0

d_all <- d %>%
  group_by(AlzID.APOE24, Study.APOE24, Data_Source, N15816_APOE24) %>%
  mutate(LabIDList.APOE24 = paste(unique(LabID.APOE24), collapse = ",")) %>%
  ungroup() %>% 
  # select(-LabID.APOE24) %>%
  distinct() %>%
  as.data.frame() 

d_all$LabIDList.APOE24 <- ifelse(d_all$AlzID.APOE24 == "-11", d_all$LabID.APOE24, d_all$LabIDList.APOE24) 

d_all <- d_all %>% 
  select(-LabID.APOE24) %>%
  distinct() %>%
  as.data.frame() 

dim(d_all)
# [1] 16007    18
length(d_all$scanName[duplicated(d_all$scanName)])

table(d_all$GWAS)
#   N Not Available            Y 
# 776          4620        10611

table(d_all$N15816_APOE24)
#   N     Y 
# 282 15725 

openxlsx::write.xlsx(d_all, "./All_Samples_from_Harmonzation_and_APOE24Table_2024_05_26.xlsx")

###--------------------------------------------------------------------------###
### output lists as required
###--------------------------------------------------------------------------###

### list 1 (N=219): Missing call rate > 5%

missing_call_rate_ex <- d_all[!is.na(d_all$GWAS_QC_Failed_Notes) & d_all$GWAS_QC_Failed_Notes == "Missing call rate > 5%",]
dim(missing_call_rate_ex)
# [1] 219  18
openxlsx::write.xlsx(missing_call_rate_ex, "./GWAS_Harmonizaton_Samples_excluded_due_to_Missing call rate > 5%_2024_05_26.xlsx")

### list 2 (N=63): Sex Mismatch

sex_mismatch_ex <- d_all[!is.na(d_all$GWAS_QC_Failed_Notes) & d_all$GWAS_QC_Failed_Notes == "mislabeled sex",]
dim(sex_mismatch_ex)
# [1] 63 18
openxlsx::write.xlsx(sex_mismatch_ex, "./GWAS_Harmonizaton_Samples_excluded_due_to_mislabeled sex_2024_05_26.xlsx")

### list 3 (N=39+1340 after checking BENIN and SLVDS): Race Discrepancies

race_apoe_discrepancies <- openxlsx::read.xlsx("./AllCohorts_AllChr_Harmonized_Race_APOE_discrepancies.xlsx")
dim(race_apoe_discrepancies)
# [1] 10611    16
colnames(race_apoe_discrepancies) <- c("scanName", "19:45411941:T:C[b37]_C", "19:45426792:G:A[b37]_A", "APOE.Imputed", "AlzID.Harmonization", "LabID.Harmonization", "StudyID", "Batch", "Sex", "Race", "APOE.Lab", "age_use_recode", "scanName2", "Case_Control_Status", "apoe_alert", "race_alert")

race_discrepancies <- race_apoe_discrepancies[!is.na(race_apoe_discrepancies$race_alert), ]
dim(race_discrepancies)
# [1] 1379   16
table(race_discrepancies$Race)
# Asian   Black Unknown   White 
#     6       8    1340      25

race_discrepancies$Study <- as.character(NA)
for (seq in 1:nrow(race_discrepancies)){
  scanName <- race_discrepancies[seq,]$scanName
  
  study.Harmonization <- ifelse(!is.na(unique(d_all[d_all$scanName %in% race_discrepancies[seq,]$scanName, ]$Study.Harmonization)), unique(d_all[d_all$scanName %in% race_discrepancies[seq,]$scanName, ]$Study.Harmonization), NA)
  study.APOE24 <- ifelse(!is.na(unique(d_all[d_all$scanName %in% race_discrepancies[seq,]$scanName, ]$Study.APOE24)), unique(d_all[d_all$scanName %in% race_discrepancies[seq,]$scanName, ]$Study.APOE24), NA)
  race_discrepancies[seq,]$Study <- ifelse(
    !is.na(study.Harmonization) & !is.na(study.APOE24) & study.Harmonization == study.APOE24,
    study.Harmonization,
    ifelse(
      is.na(study.Harmonization),
      study.APOE24,
      ifelse(
        is.na(study.APOE24),
        study.Harmonization,
        paste0(study.Harmonization, ",", study.APOE24)
      )
    )
  )
  
  race_discrepancies[seq,]$LabID.Harmonization <- ifelse(!is.na(race_discrepancies[seq,]$LabID.Harmonization), race_discrepancies[seq,]$LabID.Harmonization, unique(d_all[d_all$scanName %in% race_discrepancies[seq,]$scanName, ]$LabID.Harmonization))
}

race_discrepancies$Study[race_discrepancies$Study %in% c("ADRC,ADRC,HCP","HCP,ADRC,HCP")] <- "ADRC,HCP"
race_discrepancies$Study[race_discrepancies$Study == "SARP,BACNE,SARP"] <- "BACNE,SARP"

table(race_discrepancies$Study)
# ADRC        ADRC,HCP           BACNE      BACNE,SARP           BENIN             GEM             HCP 
#   25               2               1               1             748               4               2 
# Heartscore          IGNITE          MOVIES           MYHAT SARP,ADRC,BACNE         SENIORS           SLVDS 
#          2              10               2               6               1              18             556 
# WHITE 
#     1
sum(is.na(race_discrepancies$Study))
# 0
table(race_discrepancies$Race)
# Asian   Black Unknown   White 
#     6       8    1340      25
race_discrepancies1 <- race_discrepancies[,c("AlzID.Harmonization", "scanName", "LabID.Harmonization", "StudyID", "Study", "Batch", "Race",
                                           "Case_Control_Status", "race_alert")]
table(race_discrepancies$race_alert)
# no Race in database self-reported Asian clustered in White self-reported Black clustered in White 
#                1340                                      6                                      8 
# self-reported White clustered in Black 
#                                     25

for (seq in 1:nrow(race_discrepancies1)){
  race_discrepancies1[seq,]$Race <- ifelse(race_discrepancies1[seq,]$Race == "Unknown" & race_discrepancies1[seq,]$Study == "BENIN", "Black", 
                                     ifelse(race_discrepancies1[seq,]$Race == "Unknown" & race_discrepancies1[seq,]$Study == "SLVDS", "White", race_discrepancies1[seq,]$Race))
}

table(race_discrepancies1$Race)
# Asian   Black Unknown   White 
#     6     756      36     581

race_discrepancies1$race_alert1 <- ifelse(race_discrepancies1$race_alert %in% c("self-reported White clustered in Black", "self-reported Asian clustered in White", "self-reported Black clustered in White"), race_discrepancies1$race_alert, 
                                          ifelse(race_discrepancies1$Race == "Unknown", "no Race in database", NA))

table(race_discrepancies1$race_alert1)
# no Race in database self-reported Asian clustered in White self-reported Black clustered in White 
#                  36                                      6                                      8 
# self-reported White clustered in Black 
#                                     25

race_disc_output <- race_discrepancies1[!is.na(race_discrepancies1$race_alert1), colnames(race_discrepancies1) != "race_alert"] %>%
  rename(race_alert = race_alert1)
dim(race_disc_output)
# [1] 75 9
openxlsx::write.xlsx(race_disc_output[order(race_disc_output$race_alert),], "./GWAS_Harmonizaton_Samples_with_Self-Reported_Race_discrepancies_compared_to_PCA_2024_05_26.xlsx")

### list 4 (N = 1035+329): APOE discrepancies

apoe_discrepancies <- race_apoe_discrepancies[!is.na(race_apoe_discrepancies$apoe_alert), ]
dim(apoe_discrepancies)
# [1] 1364   16

apoe_discrepancies$Study <- as.character(NA)
for (seq in 1:nrow(apoe_discrepancies)){
  scanName <- apoe_discrepancies[seq,]$scanName
  
  study.Harmonization <- ifelse(!is.na(unique(d_all[d_all$scanName %in% apoe_discrepancies[seq,]$scanName, ]$Study.Harmonization)), unique(d_all[d_all$scanName %in% apoe_discrepancies[seq,]$scanName, ]$Study.Harmonization), NA)
  study.APOE24 <- ifelse(!is.na(unique(d_all[d_all$scanName %in% apoe_discrepancies[seq,]$scanName, ]$Study.APOE24)), unique(d_all[d_all$scanName %in% apoe_discrepancies[seq,]$scanName, ]$Study.APOE24), NA)
  apoe_discrepancies[seq,]$Study <- ifelse(
    !is.na(study.Harmonization) & !is.na(study.APOE24) & study.Harmonization == study.APOE24,
    study.Harmonization,
    ifelse(
      is.na(study.Harmonization),
      study.APOE24,
      ifelse(
        is.na(study.APOE24),
        study.Harmonization,
        paste0(study.Harmonization, ",", study.APOE24)
      )
    )
  )
  
  apoe_discrepancies[seq,]$LabID <- ifelse(!is.na(apoe_discrepancies[seq,]$LabID), apoe_discrepancies[seq,]$LabID, unique(d_all[d_all$scanName %in% apoe_discrepancies[seq,]$scanName, ]$LabID.Harmonization))
}

apoe_discrepancies$Study[apoe_discrepancies$Study %in% c("ADRC,ADRC,BACNE","BACNE,ADRC,BACNE")] <- "ADRC,BACNE"
apoe_discrepancies$Study[apoe_discrepancies$Study %in% c("ADRC,ADRC,HCP","HCP,ADRC,HCP")] <- "ADRC,HCP"
apoe_discrepancies$Study[apoe_discrepancies$Study == "BACNE,BACNE,SARP,WHITE"] <- "BACNE,SARP,WHITE"
apoe_discrepancies$Study[apoe_discrepancies$Study == "BACNE,BACNE,WHITE"] <- "BACNE,WHITE"
apoe_discrepancies$Study[apoe_discrepancies$Study == "BACNE,BACNE,SARP"] <- "BACNE,SARP"
apoe_discrepancies$Study[apoe_discrepancies$Study == "WHITE,SARP,WHITE"] <- "WHITE,SARP"
apoe_discrepancies$Study[apoe_discrepancies$Study == "GEM,GEM,BACNE,SARP"] <- "GEM,BACNE,SARP"
apoe_discrepancies$Study[apoe_discrepancies$Study == "GEM,GEM,MYHAT"] <- "GEM,MYHAT"
apoe_discrepancies$Study[apoe_discrepancies$Study == "MYHAT,MYHAT_PiB,MYHAT"] <- "MYHAT,MYHAT_PiB"

table(apoe_discrepancies$Study)
# ADRC       ADRC,BACNE         ADRC,HCP            BACNE       BACNE,SARP BACNE,SARP,WHITE      BACNE,WHITE 
#  364                2                7               23                3                1                3 
# BENIN              GEM     GEM,ADRC,GEM   GEM,BACNE,SARP        GEM,MYHAT              HCP       Heartscore 
#    47              319                1                1                3               10               15 
# IGNITE             LLDP           MOVIES            MYHAT  MYHAT,MYHAT_PiB  SARP,ADRC,BACNE          SENIORS 
#     48                2               73              364               10                1                3 
# SLVDS            WHITE       WHITE,SARP 
#    53                9                2
sum(is.na(apoe_discrepancies$Study))
# [1] 0

apoe_disc_output <- apoe_discrepancies[,c("AlzID.Harmonization", "scanName", "LabID.Harmonization", "StudyID", "Study", "Batch",
                                          "19:45411941:T:C[b37]_C", "19:45426792:G:A[b37]_A", "APOE.Imputed", "APOE.Lab", "Case_Control_Status", "apoe_alert")] %>%
  rename(`apoe_alert(Imputed|Lab)` = apoe_alert) 

for(seq in 1:nrow(apoe_disc_output)){
  
  scanName <- apoe_disc_output[seq,]$scanName
  apoe <- unique(d_all[d_all$scanName %in% scanName, ]$APOE)
  
  apoe_disc_output[seq,]$APOE.Lab <- apoe
}

for(seq in 1:nrow(apoe_disc_output)){
  apoe_disc_output[seq,]$`apoe_alert(Imputed|Lab)` <- as.character(ifelse(is.na(apoe_disc_output[seq,]$APOE.Lab), "no APOE in database",
                                                                          ifelse(apoe_disc_output[seq,]$APOE.Lab == apoe_disc_output[seq,]$APOE.Imputed, NA, paste0(apoe_disc_output[seq,]$APOE.Imputed,"|",apoe_disc_output[seq,]$APOE.Lab))))
}

apoe_disc_output1 <- apoe_disc_output[!is.na(apoe_disc_output$`apoe_alert(Imputed|Lab)`),]
sum(is.na(apoe_disc_output1$APOE.Lab))
# 313
# table(apoe_disc_output1$`apoe_alert(Imputed|Lab)`)
# tmp <- apoe_disc_output1[apoe_disc_output1$`apoe_alert(Imputed|Lab)` != "no APOE in database" & !is.na(apoe_disc_output1$`apoe_alert(Imputed|Lab)`),]
# table(tmp$`apoe_alert(Imputed|Lab)`, tmp$Batch)

openxlsx::write.xlsx(apoe_disc_output1[order(apoe_disc_output1$apoe_alert),], "./GWAS_Harmonizaton_Samples_with_APOE_discrepancies_2024_05_26.xlsx")

### list 5: differences between the APOE24 table and the GWAS list at the beginning

length(unique(d_all[d_all$GWAS == "Not Available", ]$AlzID.APOE24))
# 4577

notGWAS <- d_all[d_all$GWAS == "Not Available", ]
notGWAS$AlzID.APOE24[duplicated(notGWAS$AlzID.APOE24)]
table(notGWAS$Study.APOE24)

# ADRC      ADRC,BACNE        AgeWise4            ALPS           BACNE      BACNE,SARP           BENIN 
# 1816               1              51              31               8               1              46 
# DEKOSKY           FLOYD             GEM             HCP      Heartscore          IGNITE            LLDP 
#       5              17             127               9               4               2             522 
# MOVIES         MSBrain           MYHAT MYHAT_PiB,MYHAT             RAW            SARP           SCOOP 
#    131             262             840               1             117             163              69 
# SENIORS           SLVDS           SPOON           WHITE         WINDOWS 
#      98             199              20              34              46 

openxlsx::write.xlsx(notGWAS, "./Samples_not_GWASed_in_APOE24_Table_2024_05_26.xlsx")

length(unique(d_all[d_all$GWAS != "Not Available" & d_all$N15816_APOE24 == "N", ]$scanName))
# 281
length(unique(d_all[d_all$GWAS != "Not Available" & d_all$N15816_APOE24 == "N", ]$AlzID.Harmonization))
# 254
openxlsx::write.xlsx(d_all[d_all$GWAS != "Not Available" & d_all$N15816_APOE24 == "N", ], "./Samples_GWASed_not_in_APOE24_Table_2024_05_26.xlsx")

openxlsx::write.xlsx(d_all[!is.na(d_all$AlzID.APOE24) & !is.na(d_all$AlzID.Harmonization) & d_all$AlzID.APOE24 != d_all$AlzID.Harmonization,], 
                     "./Samples_Conflict_AlzIDs in Harmonization and in APOE24 table_2024_05_26.xlsx")

###--------------------------------------------------------------------------###
### combine the info with the autopsy table
###--------------------------------------------------------------------------###

d_all <- openxlsx::read.xlsx("./All_Samples_from_Harmonzation_and_APOE24Table_2024_05_26.xlsx")
### the autopsy file
autopsy <- openxlsx::read.xlsx("ADRC autopsy cases for Kamboh 12-2023.xlsx")
dim(autopsy)
# [1] 853   3

### the ADRC data sent by Frank (provided by Annie)
adrc <- openxlsx::read.xlsx("/Users/shiruyu/Desktop/tmp/ADRC.xlsx")
dim(adrc)
# [1] 4898   28

# remove white spaces
adrc <- data.frame(lapply(adrc, function(x) gsub("\\s+", "", x)), stringsAsFactors = FALSE)
adrc$PTIDShort_2 <- stringr::str_pad(adrc$PTIDShort, width = 4, pad = "0")
dim(adrc)
# [1] 4898   29

adrc$LabID <- ifelse(!is.na(adrc$LabIDBlood), adrc$LabIDBlood, adrc$LabIDBrain)

length(intersect(autopsy$`ADRC#`, adrc$PTIDShort_2))
# 777
diff1 <- setdiff(autopsy$`ADRC#`, adrc$PTIDShort_2)
# 75

autopsy1 <- merge(autopsy, adrc[,c("PTID","PTIDShort_2", "LabIDBlood", "LabIDBrain","LabID")], 
                  by.x = "ADRC#", by.y = "PTIDShort_2", all.x = TRUE) %>% distinct()
dim(autopsy1)
# [1] 853   7

autopsy1[autopsy1$`ADRC#` %in% autopsy1$`ADRC#`[duplicated(autopsy1$`ADRC#`)],]
# 4538

autopsy1$AlzID.Harmonization <- NA
autopsy1$AlzID.APOE24 <- NA
autopsy1$Study <- NA

for (seq in 1:nrow(autopsy1)){
  
  labID <- autopsy1[seq,]$LabID
  
  if(!is.na(labID)){
    tmp <- d_all[d_all$scanName %in% c(tolower(labID), toupper(labID)) | d_all$LabID.Harmonization %in% c(tolower(labID), toupper(labID)) | d_all$LabIDList.APOE24 %in% d_all$LabIDList.APOE24[grep(paste0("\\b",toupper(labID),"\\b"), d_all$LabIDList.APOE2)] | d_all$LabIDList.APOE24 %in% d_all$LabIDList.APOE24[grep(paste0("\\b",tolower(labID),"\\b"), d_all$LabIDList.APOE24, fixed = TRUE)],] 
    
    if(nrow(tmp) > 0){
      autopsy1[autopsy1$LabID %in% labID,]$AlzID.Harmonization <- ifelse(length(setdiff(tmp$AlzID.Harmonization, NA)) > 0, setdiff(tmp$AlzID.Harmonization, NA), NA)
      autopsy1[autopsy1$LabID %in% labID,]$AlzID.APOE24 <- ifelse(length(setdiff(tmp$AlzID.APOE24, NA)) > 0, setdiff(tmp$AlzID.APOE24, NA), NA)
      
      study.Harmonization <- unique(ifelse(!is.na(tmp$Study.Harmonization), tmp$Study.Harmonization, NA))
      study.APOE24 <- unique(ifelse(!is.na(tmp$Study.APOE24), tmp$Study.APOE24, NA))
      
      autopsy1[autopsy1$LabID %in% labID,]$Study <- ifelse(
        !is.na(study.Harmonization) & !is.na(study.APOE24) & study.Harmonization == study.APOE24,
        study.Harmonization,
        ifelse(
          is.na(study.Harmonization),
          study.APOE24,
          ifelse(
            is.na(study.APOE24),
            study.Harmonization,
            paste0(study.Harmonization, ",", study.APOE24)
          )
        )
      )
    } else {
      autopsy1[autopsy1$LabID %in% labID,]$AlzID.Harmonization <- NA
      autopsy1[autopsy1$LabID %in% labID,]$AlzID.APOE24 <- NA
      autopsy1[autopsy1$LabID %in% labID,]$Study <- NA
    }
      
  }

}

nrow(autopsy1[!is.na(autopsy1$AlzID.APOE24) | !is.na(autopsy1$AlzID.Harmonization),])
# 772
autopsy1$`ADRC#`[duplicated(autopsy1$`ADRC#`)]

openxlsx::write.xlsx(autopsy1, "./ADRC autopsy cases with AlzIDs assigned_2024_05_26.xlsx")

###--------------------------------------------------------------------------###
### save working space
###--------------------------------------------------------------------------###

# save.image("./tmp.RData")
# load("./tmp.RData")









