
rm(list = ls())
setwd("/ix/kfan/Ruyu/harmonization/All_Cohorts/")

library(data.table)
library(dplyr)

### Note: this script is used to check the data consistencies among GWAS Harmonization Table, APOE24 Table, and Autopsy 

###--------------------------------------------------------------------------###
### load required files
###--------------------------------------------------------------------------###

###--------------- the original file sent from Lynda related to harmonization
# 2023_12_20 version
oldRef_2023_12_20 <- openxlsx::read.xlsx("./Ref_files_2024_06_14/V_DataRequest_Ruyu_GWAS_Info_2023_12_20.xlsx")
length(unique(oldRef_2023_12_20$Scan_Name)) # 11159

oldRef_2023_12_20$scanName <- gsub(" ", "", oldRef_2023_12_20$Scan_Name)
# oldRef_2023_12_20$sex_recode1 <- ifelse(oldRef_2023_12_20$Sex_DB < 0, NA, oldRef_2023_12_20$Sex_DB)
# oldRef_2023_12_20$sex_recode <- ifelse(!is.na(oldRef_2023_12_20$sex_recode1), oldRef_2023_12_20$sex_recode1, oldRef_2023_12_20$Sex)
oldRef_2023_12_20$race_recode1 <- ifelse(!is.na(oldRef_2023_12_20$Race_DB), oldRef_2023_12_20$Race_DB, oldRef_2023_12_20$Race)
oldRef_2023_12_20$race_recode <- ifelse(oldRef_2023_12_20$race_recode1 < 0 | is.na(oldRef_2023_12_20$race_recode1), NA,
                             ifelse(oldRef_2023_12_20$race_recode1 %in% c("W", "White"), "White", 
                                    ifelse(oldRef_2023_12_20$race_recode1 %in% c("B", "Black"), "Black",
                                           ifelse(oldRef_2023_12_20$race_recode1 %in% c("A", "Asian"), "Asian","Others"))))

table(oldRef_2023_12_20$race_recode)
# Asian   Black  Others   White 
#    26     708      74    9169
sum(is.na(oldRef_2023_12_20$race_recode))
# 1410

# oldRef_2023_12_20$apoe_recode1 <- ifelse(oldRef_2023_12_20$APOE_DB < 0 | oldRef_2023_12_20$APOE_DB == 99 | is.na(oldRef_2023_12_20$APOE_DB), NA, oldRef_2023_12_20$APOE_DB)
# oldRef_2023_12_20$apoe_recode <- ifelse(!is.na(oldRef_2023_12_20$apoe_recode1), oldRef_2023_12_20$apoe_recode1, oldRef_2023_12_20$APOE)
oldRef_2023_12_20$case_control1 <- ifelse(oldRef_2023_12_20$Case_Control_DB == "-1" | is.na(oldRef_2023_12_20$APOE_DB), NA, oldRef_2023_12_20$Case_Control_DB)
oldRef_2023_12_20$case_control <-  ifelse(!is.na(oldRef_2023_12_20$case_control1), oldRef_2023_12_20$case_control1, oldRef_2023_12_20$Case_Control)

# 2024_05_01 version
oldRef_2024_05_01 <- openxlsx::read.xlsx("./Ref_files_2024_06_14/V_DataRequest_Ruyu_Harmonization_2024_05_01.xlsx")
oldRef_2024_05_01$scanName <- gsub(" ", "", oldRef_2024_05_01$Scan_name)
length(unique(oldRef_2024_05_01$scanName)) # 10603

oldRef_2024_05_01$race_recode <- ifelse(oldRef_2024_05_01$Race < 0 | is.na(oldRef_2024_05_01$Race), NA,
                                        ifelse(oldRef_2024_05_01$Race %in% c("W", "White"), "White", 
                                               ifelse(oldRef_2024_05_01$Race %in% c("B", "Black"), "Black",
                                                      ifelse(oldRef_2024_05_01$Race %in% c("A", "Asian"), "Asian","Others"))))
table(oldRef_2024_05_01$race_recode)
# Asian   Black  Others   White 
#    25    1361      73    9032
sum(is.na(oldRef_2024_05_01$race_recode))
# 112

oldRef_2024_05_01$case_control <-  ifelse(oldRef_2024_05_01$Case_Control_Current < 0, NA, oldRef_2024_05_01$Case_Control_Current)
table(oldRef_2024_05_01$case_control)
# Case Control    EOAD 
# 3172    7208       1

###--------------- the APOE24 table
# updated on 2024_06_11 

apoe24 <- openxlsx::read.xlsx("./Ref_files_2024_06_14/APOE24_Sample_Characteristics_SENIORS-updated_N15186_2024_06_11.xlsx")
# dim(apoe24)
# [1] 15186    32
apoe24_split <- apoe24 %>%
  tidyr::separate_rows(LabIDList, sep = ",") %>%
  as.data.frame()
# dim(apoe24_split)
# [1] 16204    32

table(apoe24_split$race_recode)
# Asian   Black  Others   White 
#    55    2031      41   13679 
sum(is.na(apoe24_split$race_recode))
# 398

table(apoe24_split$sex_recode)
#    F    M 
# 9111 6676
sum(is.na(apoe24_split$sex_recode))
# 417

apoe24_split$case_control <- ifelse(!is.na(apoe24_split$dx1_disease_status_recode), apoe24_split$dx1_disease_status_recode, apoe24_split$cc_recode)
table(apoe24_split$case_control)
#   AD         Control        Dementia Dementia-Non-AD            EOAD             MCI 
# 3310           11299             415             120               1            1035

###--------------- the harmonization phenotype file

all <- read.xlsx("./QC_Check/GWAS_Harmonization_N=11387_Demographic_Final.xlsx")
# dim(all)
# [1] 11387     9
table(all$Sex)
#    F    M 
# 6389 4998
table(all$Race)
# Asian   Black  Others Unknown   White 
#    27    1497      74      88    9701

###--------------------------------------------------------------------------###
### combine the info with the references (add phenotype) using the APOE24 table
###--------------------------------------------------------------------------###

### check LabID unique identifiers
# unique(gsub("[0-9]","",apoe24_split$LabIDList))
# unique(gsub("[0-9]","",all$V1))
table(apoe24_split$StudyList)
table(apoe24_split$studylist_recode)

# ref file received from Frank and Annie: Kamboh Lab Studies_2024-05-15.xlsx
# ADRC: A/B, -MT, -CB, -CG, -AH, -SE, UVA
# BACNE/SARP: BACN/S
# HCP: HC
# MYHAT_PiB/MYHAT: G2PiB/G2-/G (6 digits)
# Heartscore: cPiB
# AgeWise4: AW
# White: W
# RAW: RW
# IGNITE: IGN
# ALPS: ALPS (6 digits)
# FLOYD: F (6 digits)
# MSBrain: MS
# SCOOP: PCD
# WINDOWS: WBS
# SENIORS: SR
# DEKOSKY: DEKOSKY
# MOVIES: M (4 digits) : in Annie's file, lower case m is used
# LLDP: L (6 digits)
# GEM: 7 digits

# head(all)

all$LabID <- ifelse(grepl("_R", all$V1), sub("^(.*?)_.*$", "\\1", all$V1), 
                    ifelse(gsub("[0-9]","",all$V1) == "A_A", sub("^(.*?)_.*$", "\\1", all$V1), 
                           ifelse(gsub("[0-9]","",all$V1) == "B_B", sub("^(.*?)_.*$", "\\1", all$V1), all$V1)))
all$LabID <- ifelse(grepl("R$|-R$", all$LabID), gsub("R|-R", "", all$LabID), 
                           ifelse(grepl("_Kamboh_C_", all$LabID), sub(".*?(Kamboh_C_.*)", "\\1", all$LabID), all$LabID))

unique(gsub("[0-9]","",all$LabID))
all$LabID <- sub("^a", "A", all$LabID)
all$LabID <- sub("^b", "B", all$LabID)
all$LabID <- sub("^g", "G", all$LabID)
all$LabID <- sub("^s", "S", all$LabID)
all$LabID <- sub("^m", "M", all$LabID)
all$LabID <- sub("^UVa", "UVA", all$LabID)
all$LabID <- sub("^cPib", "cPiB", all$LabID)

setdiff(unique(gsub("[0-9]","",all$LabID)), unique(gsub("[0-9]","",apoe24_split$LabIDList)))
# [1] "Kamboh_C_"       "AP"        "MP"      
# weird samples: A0370P, M110P

###--------------------------------------------------------------------------###
### add the LabID for those with Kamboh_C_ prefix
###--------------------------------------------------------------------------###

# load reference files

kamboh_c_files <- list.files("./Ref_files_2024_06_14/Kamboh_C_files")
kamboh_c_list <- list()

for (file in kamboh_c_files){
  if(grepl("csv", file)){
    tmp <- read.csv(paste0("./Ref_files_2024_06_14/Kamboh_C_files/", file))
    kamboh_c_list[[file]] <- tmp[,c("labid","gwas","gwasID", "apoe", "sex", "race")]
  } else {
    tmp <- openxlsx::read.xlsx(paste0("./Ref_files_2024_06_14/Kamboh_C_files/", file))
    colnames(tmp)[colnames(tmp) == "gwasId"] <- "gwasID"
    kamboh_c_list[[file]] <- tmp[,c("labid","gwas","gwasID", "apoe", "sex", "race")]
  }
}

kamboh_c_ref <- do.call(rbind, kamboh_c_list) %>% distinct()
dim(kamboh_c_ref)
# [1] 6549    6
kamboh_c_ref$LabID <- toupper(kamboh_c_ref$labid)

# table(kamboh_c_ref$sex)
# table(kamboh_c_ref$race)
# table(kamboh_c_ref$apoe)

# replace NULL with NA
kamboh_c_ref <- as.data.frame(lapply(kamboh_c_ref, function(x) {
  x[x == "NULL"] <- NA
  return(x)
}))

# recode race
kamboh_c_ref$race_recode <- ifelse(kamboh_c_ref$race == "W", "White", 
                                   ifelse(kamboh_c_ref$race == "B", "Black",
                                          ifelse(kamboh_c_ref$race == "A", "Asian", 
                                                 ifelse(kamboh_c_ref$race == "O", "Others", NA))))
# table(kamboh_c_ref$race_recode)
# Asian  Black Others  White 
#.   12    454     10   5856
sum(is.na(kamboh_c_ref$race_recode))
# 217

# replace LabID for those with kamboh_c_ prefix in the database
# set a counter
i <- 0 
for (seq in 1:nrow(all)) {
  if (grepl("Kamboh_C_", all$LabID[seq])) {
    kamboh_id <- all$LabID[seq]
    tmp <- kamboh_c_ref[kamboh_c_ref$gwasID %in% kamboh_id, ]
    if (nrow(tmp) > 0) {  # Check if tmp has any matching rows
      kamboh_labid <- tmp$LabID
      
      # replace with the ref LabID
      all[all$LabID == kamboh_id, "LabID"] <- kamboh_labid
      # count appearence
      i <- i + 1
    }
  }
}

setdiff(unique(gsub("[0-9]","",all$LabID)), unique(gsub("[0-9]","",apoe24_split$LabIDList)))
# [1] "AP" "MP"

### merge the Harmzonization with apoe24 list by LabID
all$Description <- "GWAS Harmonization"

length(intersect(all$LabID, apoe24_split$LabIDList))
# 9294

all$Race2 <- ifelse(all$Race == "Unknown", NA, all$Race)
all$Description <- "GWAS Harmonization"
apoe24_split$Description <- "APOE24 Table"

d <- merge(all, apoe24_split[,c("AlzID", "LabIDList", "StudyList", "Description")], 
           by.x = "LabID", by.y = "LabIDList", all = TRUE) %>% distinct()
names(d)[names(d) == "V1"] <- "scanName"
# d$scanName[duplicated(d$scanName) & !is.na(d$scanName)]
d$AlzID1 <- ifelse(!is.na(d$AlzID.x) & !is.na(d$AlzID.y) & d$AlzID.x == d$AlzID.y, d$AlzID.x, 
                  ifelse(!is.na(d$AlzID.y), d$AlzID.y, d$AlzID.x))

d$AlzID <- ifelse(d$AlzID.x %in% "-12", "-12", d$AlzID1)
sum(is.na(d$AlzID))
# 0 
nrow(d[d$AlzID %in% "-11",])
# 128
nrow(d[d$AlzID %in% "-12",])
# 71
# View(d[!is.na(d$AlzID.x) & !is.na(d$AlzID) & d$AlzID.x != d$AlzID, ])
# View(apoe24_split)
# dim(d)
# [1] 18297    17
sum(is.na(d$Batch))
# 6910

d[d$scanName %in% "B0324_B0337", "LabID"] <- "B0337"
d[d$scanName %in% "B0324_B0337", "AlzID"] <- "ALZ105858"

d[d$scanName %in% "B0337_B0324", "LabID"] <- "B0324"
d[d$scanName %in% "B0337_B0324", "AlzID"] <- "ALZ105853"

d$Description <- ifelse(!is.na(d$Description.x) & !is.na(d$Description.y), "GWAS Harmonization/APOE24 Table", 
                        ifelse(!is.na(d$Description.x), d$Description.x, d$Description.y))

dd <- d[,c("AlzID","LabID", "scanName","Batch", "Removal.Reasons", "Additional_Notes",
           "GWAS_QC", "Sex", "Race2", "Study", "StudyList", "Description")] %>%
  distinct()
# dim(dd)
# [1] 18297    12

colnames(dd) <- c("AlzID","LabID", "scanName","Batch", "Removal.Reasons", "Additional_Notes",
                  "GWAS_QC", "Sex", "Race", "Study", "StudyList", "Description")

nrow(dd[dd$AlzID %in% "-11",])
# 128
nrow(dd[dd$AlzID %in% "-12",])
# 71

###--------------------------------------------------------------------------###
### fill in Sex, Race, APOE, and Case_Control Status
###--------------------------------------------------------------------------###

# Note: APOE information must all come from the final APOE24 Table sent by Lynda
#       Case_Control_Status: treat the APOE24 table as the top priority - but if the sample has been identified as "Case" in any of the old references, then assign him as "Cases"
#       Race: treat the APOE24 table as the top priority - fill in missing with older refereces
#       AlzID sticks to the AlzID_DB in the old references

dd$APOE <- as.character(NA_character_)
dd$Sex_2023_12_20 <- as.character(NA_character_)
dd$Sex_2024_05_01 <- as.character(NA_character_)
dd$Race_2023_12_20 <- as.character(NA_character_)
dd$Race_2024_05_01 <- as.character(NA_character_)
dd$CC_2023_12_20 <- as.character(NA_character_)
dd$CC_2024_05_01 <- as.character(NA_character_)
dd$CC_APOE24 <- as.character(NA_character_)

for (seq in 1:nrow(dd)){
  
  # step 1: using the APOE24 to fill in information of interest
  # sex
  if (is.na(dd[seq, "Sex"]) & dd[seq, "AlzID"] != "-11"){
    matching_row <- apoe24_split[apoe24_split$AlzID %in% dd[seq, "AlzID"],]
    if(nrow(matching_row) > 0){
      dd[seq, "Sex"] <- unique(matching_row$sex_recode)
    }
  }
  if (is.na(dd[seq, "Sex"])){
    matching_row <- apoe24_split[apoe24_split$LabIDList %in% dd[seq, "LabID"],]
    if(nrow(matching_row) > 0){
      dd[seq, "Sex"] <- unique(matching_row$sex_recode)
    }
  }
  # race
  if (is.na(dd[seq, "Race"]) & dd[seq, "AlzID"] != "-11"){
    matching_row <- apoe24_split[apoe24_split$AlzID %in% dd[seq, "AlzID"],]
    if(nrow(matching_row) > 0){
      dd[seq, "Race"] <- unique(matching_row$race_recode)
    }
  }
  if (is.na(dd[seq, "Race"])){
    matching_row <- apoe24_split[apoe24_split$LabIDList %in% dd[seq, "LabID"],]
    if(nrow(matching_row) > 0){
      dd[seq, "Race"] <- unique(matching_row$race_recode)
    }
  }
  # APOE
  if (is.na(dd[seq, "APOE"]) & !(dd[seq, "AlzID"] %in% c("-11","-12"))){
    matching_row <- apoe24_split[apoe24_split$AlzID %in% dd[seq, "AlzID"],]
    if(nrow(matching_row) > 0){
      dd[seq, "APOE"] <- unique(matching_row$apoe_recode)
    }
  }
  if (is.na(dd[seq, "APOE"])){
    matching_row <- apoe24_split[apoe24_split$LabIDList %in% dd[seq, "LabID"],]
    if(nrow(matching_row) > 0){
      dd[seq, "APOE"] <- unique(matching_row$apoe_recode)
    }
  }
  # Case_Control_Status
  if (is.na(dd[seq, "CC_APOE24"]) & dd[seq, "AlzID"] != "-11"){
    matching_row <- apoe24_split[apoe24_split$AlzID %in% dd[seq, "AlzID"],]
    if(nrow(matching_row) > 0){
      dd[seq, "CC_APOE24"] <- unique(matching_row$case_control)
    }
  }
  if (is.na(dd[seq, "CC_APOE24"])){
    matching_row <- apoe24_split[apoe24_split$LabIDList %in% dd[seq, "LabID"],]
    if(nrow(matching_row) > 0){
      dd[seq, "CC_APOE24"] <- unique(matching_row$case_control)
    }
  }
  
  # step 2: using the ref files to fill in information of interest
  if (is.na(dd[seq, "Sex_2023_12_20"])){
    matching_row <- oldRef_2023_12_20[oldRef_2023_12_20$scanName %in% dd[seq, "scanName"],]
    if(nrow(matching_row) > 0){
      dd[seq, "Sex_2023_12_20"] <- ifelse(length(unique(matching_row$Sex_DB)) > 0, 
                                          unique(matching_row$Sex_DB)[!is.na(unique(matching_row$Sex_DB))], unique(matching_row$Sex_DB))
    }
  }
  if (is.na(dd[seq, "Sex_2024_05_01"])){
    matching_row <- oldRef_2024_05_01[oldRef_2024_05_01$scanName %in% dd[seq, "scanName"],]
    if(nrow(matching_row) > 0){
      dd[seq, "Sex_2024_05_01"] <- ifelse(length(unique(matching_row$Sex)) > 0, 
                                          unique(matching_row$Sex)[!is.na(unique(matching_row$Sex))], unique(matching_row$Sex))
    }
  }
  
  if (is.na(dd[seq, "Race_2023_12_20"])){
    matching_row <- oldRef_2023_12_20[oldRef_2023_12_20$scanName %in% dd[seq, "scanName"],]
    if(nrow(matching_row) > 0){
      dd[seq, "Race_2023_12_20"] <- ifelse(length(unique(matching_row$race_recode)) > 0, 
                                           unique(matching_row$race_recode)[!is.na(unique(matching_row$race_recode))], unique(matching_row$race_recode))
    }
  }
  if (is.na(dd[seq, "Race_2024_05_01"])){
    matching_row <- oldRef_2024_05_01[oldRef_2024_05_01$scanName %in% dd[seq, "scanName"],]
    if(nrow(matching_row) > 0){
      dd[seq, "Race_2024_05_01"] <- ifelse(length(unique(matching_row$race_recode)) > 0, 
                                           unique(matching_row$race_recode)[!is.na(unique(matching_row$race_recode))], unique(matching_row$race_recode))
    }
  }
  
  if (is.na(dd[seq, "CC_2023_12_20"])){
    matching_row <- oldRef_2023_12_20[oldRef_2023_12_20$scanName %in% dd[seq, "scanName"],]
    if(nrow(matching_row) > 0){
      dd[seq, "CC_2023_12_20"] <- ifelse(length(unique(matching_row$Case_Control_DB)) > 0, 
                                         unique(matching_row$Case_Control_DB)[!is.na(unique(matching_row$Case_Control_DB))], unique(matching_row$Case_Control_DB))
    }
  }
  if (is.na(dd[seq, "CC_2024_05_01"])){
    matching_row <- oldRef_2024_05_01[oldRef_2024_05_01$scanName %in% dd[seq, "scanName"],]
    if(nrow(matching_row) > 0){
      dd[seq, "CC_2024_05_01"] <- ifelse(length(unique(matching_row$Case_Control_Current)) > 0, 
                                         unique(matching_row$Case_Control_Current)[!is.na(unique(matching_row$Case_Control_Current))], unique(matching_row$Case_Control_Current))
    }
  }
}

# recode sex, race and case_control status
dd$sex_recode <- ifelse(!is.na(dd$Sex), dd$Sex, 
                        ifelse(!is.na(dd$Sex_2023_12_20), dd$Sex_2023_12_20, dd$Sex_2024_05_01))
table(dd$sex_recode)
#     F     M 
# 10384  7549 
sum(is.na(dd$sex_recode))
# 364
# table(dd[!is.na(dd$Batch), ]$sex_recode)
#    F    M 
# 6389 4998

dd$race_recode <- ifelse(!is.na(dd$Race), dd$Race, 
                        ifelse(!is.na(dd$Race_2023_12_20), dd$Race_2023_12_20, dd$Race_2024_05_01))
table(dd$race_recode)
# Asian  Black Others  White 
#    56   2191     90  15559
sum(is.na(dd$race_recode))
# 401

dd$cc_recode <- as.character(NA_character_)
for (seq in 1:nrow(dd)){
  dd[seq, "cc_recode"] <- ifelse(grepl("AD|Case|Dementia|MCI", dd[seq, "CC_APOE24"]) | grepl("AD|Case|Dementia|MCI", dd[seq, "CC_2023_12_20"]) | grepl("AD|Case|Dementia|MCI", dd[seq, "CC_2024_05_01"]), "Case", 
                                 ifelse(!is.na(dd[seq, "CC_APOE24"]), dd[seq, "CC_APOE24"], 
                                        ifelse(!is.na(dd[seq, "CC_2023_12_20"]), dd[seq, "CC_2023_12_20"], dd[seq, "CC_2024_05_01"])))
}

table(dd$cc_recode)
dd$Case_Control_Status <- ifelse(dd$cc_recode < 0 | is.na(dd$cc_recode), NA, 
                                 dd$cc_recode)
dd[dd$AlzID == "ALZ115938", "Case_Control_Status"] <- "Case"
table(dd$Case_Control_Status)
# Case Control 
# 5271   13006
sum(is.na(dd$Case_Control_Status))
# 20

table(dd$APOE)
# 22    23    24    33    34    44 
# 98  1921   410 10299  4598   688 
sum(is.na(dd$APOE))
# 283

# recode studylist
dd$studylist_recode <- as.character(NA_character_)
for (seq in 1:nrow(dd)){
  dd[seq, "studylist_recode"] <- ifelse(!is.na(dd[seq, "Study"]) & !is.na(dd[seq, "StudyList"]) & grepl(dd[seq, "Study"], dd[seq, "StudyList"]), dd[seq, "StudyList"], 
                                        ifelse(!is.na(dd[seq, "StudyList"]), dd[seq, "StudyList"], dd[seq, "Study"]))
}

table(dd$studylist_recode)
sum(is.na(dd$studylist_recode))
# 0

###--------------------------------------------------------------------------###
### group by LabID and AlzID
###--------------------------------------------------------------------------###

d <- dd[,c("AlzID","LabID","scanName","Batch","Removal.Reasons","Additional_Notes","GWAS_QC",
           "sex_recode","race_recode","APOE","Case_Control_Status","studylist_recode","Description")]
colnames(d) <- c("AlzID","LabID","scanName","Batch","Removal.Reasons","Additional_Notes","GWAS_QC",
                 "Sex","Race","APOE","Case_Control_Status","StudyList","Description")
dim(d)
# 18297    13

d1 <- d %>%
  group_by(LabID) %>%
  mutate(
    Sex = if_else(is.na(LabID), Sex, coalesce(Sex, dplyr::first(Sex[!is.na(Sex)]))),
    Race = if_else(is.na(LabID), Race, coalesce(Race, dplyr::first(Race[!is.na(Race)]))),
    APOE = if_else(is.na(LabID), APOE, coalesce(APOE, dplyr::first(APOE[!is.na(APOE)]))),
    Case_Control_Status = if_else(is.na(LabID), Case_Control_Status, coalesce(Case_Control_Status, dplyr::first(Case_Control_Status[!is.na(Case_Control_Status)]))),
    ) %>%
  ungroup()

d1_noAlzID <- d1[d1$AlzID < 0, ]
d1_AlzID <- d1[d1$AlzID > 0,]

dd_AlzID <- d1_AlzID %>%
  group_by(AlzID) %>%
  mutate(
    Sex = if_else(is.na(AlzID), Sex, coalesce(Sex, dplyr::first(Sex[!is.na(Sex)]))),
    Race = if_else(is.na(AlzID), Race, coalesce(Race, dplyr::first(Race[!is.na(Race)]))),
    APOE = if_else(is.na(AlzID), APOE, coalesce(APOE, dplyr::first(APOE[!is.na(APOE)]))),
    Case_Control_Status = if_else(is.na(AlzID), Case_Control_Status, coalesce(Case_Control_Status, dplyr::first(Case_Control_Status[!is.na(Case_Control_Status)]))),
  ) %>%
  ungroup()

dd_grouped <- dd_AlzID %>%
  group_by(AlzID) %>%
  mutate(LabIDList = paste(unique(LabID), collapse = ",")) %>%
  ungroup() %>%
  distinct() %>%
  as.data.frame() 

d1_noAlzID$LabIDList <- d1_noAlzID$LabID
dd <- rbind(d1_noAlzID, dd_grouped) %>% as.data.frame()

### get unique LabIDList values to count # overlaps in Harmonization and APOE24
d_byLabIDList <- data.frame(LabIDList = unique(dd$LabIDList))

d_byLabIDList$AlzID <- as.character(NA_character_)
d_byLabIDList$Description <- as.character(NA_character_)
for (seq in 1:nrow(d_byLabIDList)){
  
  ids <- d_byLabIDList[seq, "LabIDList"]
  
  tmp_grouped <- dd[dd$LabIDList %in% ids, ]
  if (length(grep("-12", unique(tmp_grouped$AlzID))) > 0) {
    d_byLabIDList[seq,"AlzID"] <- "-12"
  } else {
    d_byLabIDList[seq,"AlzID"] <- unique(tmp_grouped$AlzID)
  }
  
  d_byLabIDList[seq,"Description"] <- unique(ifelse(length(unique(tmp_grouped$Description)) > 1, "GWAS Harmonization/APOE24 Table", unique(tmp_grouped$Description)))
}

dim(d_byLabIDList)
# [1] 15502     3
table(d_byLabIDList$Description)
# APOE24 Table 
#         4644
# GWAS Harmonization 
#                281 
# GWAS Harmonization/APOE24 Table 
#                           10577

final <- merge(d_byLabIDList, dd[,c("LabID","scanName","Batch","Removal.Reasons",
                                          "Additional_Notes","GWAS_QC","Sex","Race","APOE","Case_Control_Status",
                                          "StudyList", "LabIDList")], 
               by = "LabIDList") %>% distinct()
dim(final)
# [1] 18297    14

SENIORS2 <- read.xlsx("./Ref_files_2024_06_14/SENIORS_Demographic_2024_06_11.xlsx")
final[final$LabIDList == "SR004", ]$Race <- SENIORS2[SENIORS2$LabIDList == "SR004", ]$race_recode

###--------------------------------------------------------------------------###
### output the table
###--------------------------------------------------------------------------###

# table(final$GWAS_QC)
### update LabIDs for SLVDS
SLVDS_final <- final[grepl("SLVDS", final$StudyList),]
SLVDS_final$updated_LabID <- as.character(NA_character_)
for (seq in 1:nrow(SLVDS_final)){
    SLVDS_final[seq, "updated_LabID"] <- ifelse(nchar(SLVDS_final[seq, "LabID"]) < 4, paste0("0", SLVDS_final[seq, "LabID"]), SLVDS_final[seq, "LabID"])
}
SLVDS_final$LabIDList <- SLVDS_final$updated_LabID
SLVDS_final2 <- SLVDS_final[, c("LabIDList","AlzID","Description",
                                "updated_LabID","scanName","Batch","Removal.Reasons","Additional_Notes","GWAS_QC","Sex","Race", 
                                "APOE","Case_Control_Status","StudyList")]
colnames(SLVDS_final2) <- c("LabIDList","AlzID","Description",
                            "LabID","scanName","Batch","Removal.Reasons","Additional_Notes","GWAS_QC","Sex","Race", 
                            "APOE","Case_Control_Status","StudyList")
final2 <- rbind(final[!grepl("SLVDS", final$StudyList),], SLVDS_final2)
dim(final2)
# [1] 18297    14
write.xlsx(final2, "./combined_Harmonization_APOE24Table_2024_06_15.xlsx")

final2$studylist_recode <- ifelse(grepl("ADRC", final2$StudyList), "ADRC", 
                              ifelse(final2$StudyList %in% c("BACNE,SARP", "BACNE,SARP,WHITE", "BACNE,WHITE", "GEM,BACNE,SARP", "GEM,BACNE,WHITE"), "BACNE", 
                                     ifelse(final2$StudyList == "MYHAT_PiB,MYHAT", "MYHAT_PiB",
                                            ifelse(final2$StudyList == "GEM,MYHAT", "MYHAT", 
                                                   ifelse(final2$StudyList == "SARP,WHITE", "SARP", final2$StudyList)))))
sum(is.na(final2$StudyList))
### extract rows in the GWAS Harmonization only
d_harmonization <- final2[!is.na(final2$GWAS_QC),] %>% distinct()

table(d_harmonization$Sex)
#    F    M 
# 6389 4998
sum(is.na(d_harmonization$Sex))
# 0

table(d_harmonization$Race)
# Asian  Black Others  White 
#    27   1512     74   9758 
sum(is.na(d_harmonization$Race))
# 16

table(d_harmonization$Case_Control_Status)
#   AD Control 
# 3723    7663
sum(is.na(d_harmonization$Case_Control_Status))
# 1

table(d_harmonization$APOE)
# 22   23   24   33   34   44 
# 54 1130  270 6083 3068  500 
sum(is.na(d_harmonization$APOE))
# 282
# d_harmonization_noAPOE <- d_harmonization[is.na(d_harmonization$APOE),]
# unique(d_harmonization_noAPOE$Description)
# table(d_harmonization$Removal.Reasons)
table(d_harmonization$Batch, d_harmonization$StudyList)
write.xlsx(d_harmonization, "./Harmonization_Phenotype_2024_06_15.xlsx")

### extract rows passed GWAS only
d_harmonization_QCed <- final2[final2$GWAS_QC == "Pass" & !is.na(final2$scanName),] %>% distinct()

table(d_harmonization_QCed$Sex)
#    F    M 
# 5877 4554
sum(is.na(d_harmonization_QCed$Sex))
# 0

table(d_harmonization_QCed$Race)
# Asian  Black Others  White 
#    25   1405     74   8911 
sum(is.na(d_harmonization_QCed$Race))
# 16

table(d_harmonization_QCed$Case_Control_Status)
#   AD Control 
# 3334    7097 
sum(is.na(d_harmonization_QCed$Case_Control_Status))
# 0

table(d_harmonization_QCed$APOE)
# 22   23   24   33   34   44 
# 49 1043  246 5607 2795  452 
sum(is.na(d_harmonization_QCed$APOE))
# 265
# d_harmonization_noAPOE <- d_harmonization_QCed[is.na(d_harmonization_QCed$APOE),]

length(unique(d_harmonization_QCed$LabID))
nrow(d_harmonization_QCed[d_harmonization_QCed$AlzID %in% c("-11","-12"),])
# 141
table(d_harmonization_QCed$Batch, d_harmonization_QCed$StudyList)
write.xlsx(d_harmonization_QCed, "./Harmonization_QC_Passed_Phenotype_2024_06_15.xlsx")

### extract rows in APOE24 only - as a double check
d_APOE24 <- final2[is.na(final2$GWAS_QC),] %>% distinct()
dim(d_APOE24)
# 6910   15
table(d_APOE24$studylist_recode)
# d_APOE24_dups <- d_APOE24[d_APOE24$LabID %in% d_APOE24$LabID[duplicated(d_APOE24$LabID)],]

###--------------------------------------------------------------------------###
### data visualization
###--------------------------------------------------------------------------###

library(eulerr)
library(venn)
library(RColorBrewer)
library(ggplot2)
library(ggprism)

# create an empty list
# length(unique(dd$LabID))
two_groups <- c("GWAS Harmonization", "APOE24 Table")
venn_prep <- euler(c("GWAS Harmonization" = nrow(d_byLabIDList[d_byLabIDList$Description == "GWAS Harmonization",]), 
                     "APOE24 Table" = nrow(d_byLabIDList[d_byLabIDList$Description == "APOE24 Table",]), 
                     "GWAS Harmonization&APOE24 Table" = nrow(d_byLabIDList[d_byLabIDList$Description == "GWAS Harmonization/APOE24 Table",])))

pdf("./combined_Harmonization_APOE24Table_vennPlot_2024_06_15.pdf", width = 10, height = 6)
plot(venn_prep,
     fills = list(fill = c("#FFABAB", "#AFCBFF"), alpha = 0.5),
     labels = list(col = "black"), 
     quantities = TRUE)
dev.off()

### update the list with the updated LabIDs
d_byLabIDList0 <- data.frame(LabIDList = unique(final2$LabIDList))

for (seq in 1:nrow(d_byLabIDList0)){
  
  ids <- d_byLabIDList0[seq, "LabIDList"]
  
  tmp_grouped <- final2[final2$LabIDList %in% ids, ]
  d_byLabIDList0[seq,"AlzID"] <- unique(tmp_grouped$AlzID)
  
  d_byLabIDList0[seq,"Description"] <- unique(ifelse(length(unique(tmp_grouped$Description)) > 1, "GWAS Harmonization/APOE24 Table", unique(tmp_grouped$Description)))
}

dim(d_byLabIDList0)
# [1] 15502     3
table(d_byLabIDList0$Description)

##write.xlsx(d_byLabIDList0[order(d_byLabIDList0$Description),], "./combined_Harmonization_APOE24Table_grouped_AlzIDs_2024_06_15.xlsx")


dd_APOE24 <- final2[final2$LabIDList %in% d_byLabIDList0[d_byLabIDList0$Description == "APOE24 Table",]$LabIDList, 
                    c("LabIDList","AlzID","Description","Sex","Race","APOE","Case_Control_Status","StudyList","studylist_recode")] %>% distinct()

write.xlsx(dd_APOE24, "./Samples_not_GWASed_in_APOE24_Table_2024_06_15.xlsx")

dd_Harmonization <- final2[final2$LabIDList %in% d_byLabIDList0[d_byLabIDList0$Description == "GWAS Harmonization",]$LabIDList, 
                    c("LabIDList","AlzID","Description","Sex","Race","APOE","Case_Control_Status","StudyList","studylist_recode")] %>% distinct()
write.xlsx(dd_Harmonization, "./Samples_GWASed_not_in_APOE24_Table_2024_06_15.xlsx")

###--------------------------------------------------------------------------###
### save working space
###--------------------------------------------------------------------------###

# save.image("./Harmonization_Pheno_2024_06_15.RData")
# load("./Harmonization_Pheno_2024_06_15.RData")









