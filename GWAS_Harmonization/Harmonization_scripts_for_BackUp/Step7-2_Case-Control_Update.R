
library(data.table)
library(dplyr)

setwd("/ix/kfan/Ruyu/harmonization/All_Cohorts")

###--------------------------------------------------------------------------###
### load data
###--------------------------------------------------------------------------###

### read the Harmonization list with AlzID
d <- fread("./All_Cohorts_2nd_Imputation_ID_Batch_with_Sex_Race.txt")

### Process the 15186 list
db <- openxlsx::read.xlsx("./APOE24_Sample_Characteristics_for_TabPrep_N15186.xlsx")

db_split <- db %>%
  tidyr::separate_rows(LabIDList, sep = ",") %>%
  as.data.frame()
# head(db_split)

# APOE
table(db_split$APOE)
db_split$apoe_recode <- ifelse(db_split$APOE < 0 | is.na(db_split$APOE), NA, db_split$APOE)
table(db_split$apoe_recode)

# Age_To_Use
summary(db_split$Age_To_Use)
db_split$age_use_recode <- ifelse(db_split$Age_To_Use < 0 | is.na(db_split$Age_To_Use), NA, db_split$Age_To_Use)
summary(db_split$age_use_recode) 

###--------------------------------------------------------------------------###
### check overlapped IDs
###--------------------------------------------------------------------------###

length(unique(d$AlzID))
length(unique(db_split$AlzID))

length(intersect(d$AlzID, db_split$AlzID))
length(setdiff(d$AlzID, db_split$AlzID))

###--------------------------------------------------------------------------###
### get APOE and age info from APOE24 reference file
###--------------------------------------------------------------------------###

d1 <- merge(d[!is.na(d$AlzID),], db_split[!is.na(db_split$AlzID),c("AlzID","apoe_recode","age_use_recode")], by = "AlzID", all.x = TRUE) %>%
  distinct()

length(unique(d1$AlzID))
# 10489
# head(d1)

d2 <- merge(d[is.na(d$AlzID),], db_split[is.na(db_split$AlzID),c("AlzID","LabIDList","apoe_recode","age_use_recode")], 
            by.x = "LabID", by.y = "LabIDList", all.x = TRUE) %>%
  distinct()
# head(d2)
d2 <- d2[,c("AlzID.x", "scanName", "LabID", "StudyID", "Batch", "sex_recode", "race_recode", "apoe_recode", "age_use_recode")]
colnames(d2) <- c("AlzID", "scanName", "LabID", "StudyID", "Batch", "sex_recode", "race_recode", "apoe_recode", "age_use_recode")

d <- rbind(d1, d2)
dim(d)

sum(is.na(d$sex_recode))
# 0
nrow(d[d$race_recode == "Unknown",])
# 1340
# sum(is.na(d$race_recode))

# perform another round of value filling, given somoe AlzIDs are not match for the same LabID, the previous merge may miss some
for (seq in 1:nrow(d)){
  ID <- d[seq,]$scanName
  
  tmp <- db_split[db_split$LabIDList %in% c(toupper(ID),tolower(ID)),]
  if(nrow(tmp) > 0){

    d[d$scanName == ID,]$race_recode <- ifelse(is.na(d[d$scanName == ID,]$race_recode), tmp$race_recode, d[d$scanName == ID,]$race_recode)
    d[d$scanName == ID,]$apoe_recode <- ifelse(is.na(d[d$scanName == ID,]$apoe_recode), tmp$apoe_recode, d[d$scanName == ID,]$apoe_recode)
    d[d$scanName == ID,]$age_use_recode <- ifelse(is.na(d[d$scanName == ID,]$age_use_recode), tmp$age_use_recode, d[d$scanName == ID,]$age_use_recode)
  }
}

d2 <- d %>%
  group_by(LabID) %>%
  mutate(
    sex_recode = if_else(is.na(LabID), sex_recode, coalesce(sex_recode, first(sex_recode[!is.na(sex_recode)]))),
    race_recode = if_else(is.na(LabID), race_recode, coalesce(race_recode, first(race_recode[!is.na(race_recode)]))),
    apoe_recode = if_else(is.na(LabID), apoe_recode, coalesce(apoe_recode, first(apoe_recode[!is.na(apoe_recode)]))),
    age_use_recode = if_else(is.na(LabID), age_use_recode, coalesce(age_use_recode, first(age_use_recode[!is.na(age_use_recode)])))
  ) %>%
  ungroup()

d3 <- d2 %>%
  group_by(AlzID) %>%
  mutate(
    sex_recode = if_else(is.na(AlzID), sex_recode, coalesce(sex_recode, first(sex_recode[!is.na(sex_recode)]))),
    race_recode = if_else(is.na(AlzID), race_recode, coalesce(race_recode, first(race_recode[!is.na(race_recode)]))),
    apoe_recode = if_else(is.na(AlzID), apoe_recode, coalesce(apoe_recode, first(apoe_recode[!is.na(apoe_recode)]))),
    age_use_recode = if_else(is.na(AlzID), age_use_recode, coalesce(age_use_recode, first(age_use_recode[!is.na(age_use_recode)])))
  ) %>%
  ungroup()

sum(is.na(d$apoe_recode))
# 320
sum(is.na(d2$apoe_recode))
# 315
sum(is.na(d3$apoe_recode))
# 315

sum(is.na(d$age_use_recode))
# 551
sum(is.na(d2$age_use_recode))
# 551
sum(is.na(d3$age_use_recode))
# 551

###--------------------------------------------------------------------------###
### load Case/Control info
###--------------------------------------------------------------------------###

cases <- fread("./pseudoGWAS_input/AllCohorts_Cases_ID_forPlink.txt")
controls <- fread("./pseudoGWAS_input/All_Cohorts_Control_Phenotypes_for_PLINK.txt")

d3$scanName2 <- paste0(paste0(d3$scanName, "_", d3$scanName))
d3$case_control <- NA
for (seq in 1:nrow(d3)){
  d3$case_control[seq] <- ifelse(length(grep(d3$scanName2[seq], cases$IID)) > 0, "Case", "Control")
}

table(d3$case_control)
# Case Control 
# 3276    7477 

# check missing in Case/Control
table(d3$case_control, is.na(d3$apoe_recode))

#         FALSE TRUE
# Case     3164  101
# Control  7257  214

table(d3$case_control, is.na(d3$age_use_recode))

#         FALSE TRUE
# Case     3127  149
# Control  7075  402

table(d3$case_control, d3$race_recode)
#         Asian Black Others Unknown White
# Case        8   180      9       0  3079
# Control    18   473     65    1340  5581

write.table(d3, "./All_Cohorts_2nd_Imputation_Phenotype.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

################################################################################
### update Case/Control
################################################################################

# 2/1 case/control

d3$case_control_recode <- ifelse(d3$case_control == "Case", 2, 
                                ifelse(d3$case_control == "Control", 1, NA))
table(d3$case_control_recode)
sum(is.na(d3$case_control_recode))

annot_fam <- fread("./AllCohorts_AllChr_IDcorrected_SexAdded_Annot.fam")
# head(annot_fam)
length(intersect(annot_fam$V1, d3$scanName))

annot_fam2 <- merge(annot_fam[,c(paste0("V",1:5))], d3[,c("scanName", "case_control_recode")], by.x = "V2", by.y = "scanName") %>%
  arrange(match(V1, annot_fam$V1))

colnames(annot_fam2) <- c(paste0("V",1:6))
# head(annot_fam2)
write.table(annot_fam2, "./AllCohorts_AllChr_IDcorrected_SexAdded_Annot.fam", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

