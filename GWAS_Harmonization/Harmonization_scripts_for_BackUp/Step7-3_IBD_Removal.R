
library(dplyr)
library(data.table)
library(openxlsx)

###--------------------------------------------------------------------------###
### check demographics of IBD pairs identified in Harmonization Project
### 
###   2024/06/10
###   get AlzID, sex and race information for each pair of samples
###   using the old reference file sent from Lynda 
###     split the file into two sheets using pi_hat >= 0.8 as a threshold
###--------------------------------------------------------------------------###

setwd("/ix/kfan/Ruyu/harmonization/All_Cohorts/")

###--------------------------------------------------------------------------###
### load IBD information
###--------------------------------------------------------------------------###

ibd <- fread("./QC_Check/AllCohorts_AllChr_IBD.genome")
dim(ibd)
# 630  14
nrow(ibd[ibd$PI_HAT >= 0.8,])
# 152

###--------------------------------------------------------------------------###
### load demographic file
###--------------------------------------------------------------------------###

all <- openxlsx::read.xlsx("./QC_Check/GWAS_Harmonization_N=11387_Demographic_Final.xlsx")

###--------------------------------------------------------------------------###
### get required demographic info
###--------------------------------------------------------------------------###

ibd1 <- merge(ibd[,c("FID1", "IID1", "FID2", "IID2", "PI_HAT")], 
              all[all$GWAS_QC == "Pass",c("V1","AlzID","Study","Batch","Sex","Race")], 
              by.x = "IID1", by.y = "V1", all.x = TRUE) %>% distinct()

colnames(ibd1) <- c("IID1","FID1","FID2","IID2","PI_HAT","AlzID1","Study1","Batch1","Sex1","Race1")

ibd2 <- merge(ibd1, 
              all[all$GWAS_QC == "Pass",c("V1","AlzID","Study","Batch","Sex","Race")], 
              by.x = "IID2", by.y = "V1", all.x = TRUE) %>% distinct()

colnames(ibd2) <- c("IID2","IID1","FID1","FID2","PI_HAT","AlzID1","Study1","Batch1","Sex1","Race1",
                    "AlzID2","Study2","Batch2","Sex2","Race2")

ibd_check <- ibd2[order(ibd2$PI_HAT, decreasing = TRUE), 
                  c("PI_HAT","FID1","IID1","AlzID1","Study1","Batch1","Sex1","Race1",
                    "FID2","IID2","AlzID2","Study2","Batch2","Sex2","Race2")]

ibd_check$Notes <- as.character(NA_character_)
for (seq in 1:nrow(ibd_check)){
  ibd_check[seq, "Notes"] <- ifelse(grepl(toupper(ibd_check[seq,"IID1"]), toupper(ibd_check[seq,"IID2"])) | grepl(toupper(ibd_check[seq,"IID2"]), toupper(ibd_check[seq,"IID1"])) & 
                                      ibd_check[seq,"AlzID1"] %in% ibd_check[seq,"AlzID2"], "Identical IID and AlzID", 
                                    ifelse(ibd_check[seq,"AlzID1"] %in%  ibd_check[seq,"AlzID2"], "Different IID but Identical AlzID", 
                                           ifelse(!(ibd_check[seq,"AlzID1"] %in% ibd_check[seq,"AlzID2"]) & ibd_check[seq,"PI_HAT"] > 0.8, "Different IID and AlzID", NA)))
}

###--------------------------------------------------------------------------###
### output file
###--------------------------------------------------------------------------###

ibd_thre1 <- ibd_check[ibd_check$PI_HAT >= 0.8,]
dim(ibd_thre1)
# [1] 152  16
table(ibd_thre1$Notes)
# Different IID and AlzID Different IID but Identical AlzID           Identical IID and AlzID 
#                      98                                21                                33  

ibd_thre2 <- ibd_check[ibd_check$PI_HAT < 0.8,]
dim(ibd_thre2)
# [1] 478  16
table(ibd_thre2$Notes)
# Different IID but Identical AlzID           Identical IID and AlzID 
#                                35                                89

##wb <- createWorkbook()

##addWorksheet(wb, "PI_HAT >= 0.8")
##writeData(wb, "PI_HAT >= 0.8", ibd_thre1)

##addWorksheet(wb, "PI_HAT < 0.8")
##writeData(wb, "PI_HAT < 0.8", ibd_thre2)

##saveWorkbook(wb, "./QC_Check/2024_06_10/IBD_Check_2024_06_10.xlsx", overwrite = TRUE)

###--------------------------------------------------------------------------###
### get removal list of IBD check output
###--------------------------------------------------------------------------###

library(readxl)

### set priority order and assign "KEEP" to samples
# set priority order
priority_order <- c("GWAS5", "GEM", "GWAS4", "GWAS3", "GWAS2", "GWAS1")

# write a function to exclude samples
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

### sheet 1
ibd_gt_0.8 <- read_excel("./QC_Check/IBD_Check_2024_06_10_Frank.xlsx", sheet = "PI_HAT >= 0.8")
# note: in the original GWAS1, two samples have flip FID and IID, based on the IBD output, trust the IID is the correct one
# ibd_gt_0.8[ibd_gt_0.8$IID1 == "B0263_B0264", "AlzID1"] # ALZ106437
# ibd_gt_0.8[ibd_gt_0.8$IID1 == "B0264_B0263", "AlzID1"] # ALZ106058

ibd_gt_0.8[ibd_gt_0.8$IID1 == "B0263_B0264", "AlzID1"] <- "ALZ106058"
ibd_gt_0.8[ibd_gt_0.8$IID1 == "B0264_B0263", "AlzID1"] <- "ALZ106437"
ibd_gt_0.8[ibd_gt_0.8$IID1 %in% c("B0263_B0264","B0264_B0263"), "Notes"] <- "Identical IID and AlzID"

# apply the function
ibd_gt_0.8_processed <- ibd_gt_0.8 %>%
  rowwise() %>%
  mutate(Keep = determine_keep(Batch1, Batch2, priority_order)) %>%
  ungroup()

# remove samples that does not make sense
ibd_gt_0.8_processed$Keep2 <- as.character(NA_character_)
ibd_gt_0.8_processed$Remove <- as.character(NA_character_)

for (seq in 1:nrow(ibd_gt_0.8_processed)){
  ibd_gt_0.8_processed[seq, "Keep2"] <- ifelse(ibd_gt_0.8_processed[seq, "Sex"] == 0, setdiff(unlist(unique(c(ibd_gt_0.8_processed[seq, "Batch1"], ibd_gt_0.8_processed[seq, "Batch2"]))), "GWAS1"), 
                                               ibd_gt_0.8_processed[seq, "Keep"])
  ibd_gt_0.8_processed[seq, "Remove"] <- ifelse(ibd_gt_0.8_processed[seq, "Sex"] == 0, "GWAS1", 
                                                setdiff(unlist(unique(c(ibd_gt_0.8_processed[seq, "Batch1"], ibd_gt_0.8_processed[seq, "Batch2"]))), ibd_gt_0.8_processed[seq, "Keep"]))
}
table(ibd_gt_0.8_processed$Keep2)
sum(is.na(ibd_gt_0.8_processed$Keep2))
# 3

table(ibd_gt_0.8_processed$Remove)
sum(is.na(ibd_gt_0.8_processed$Remove))
# 15 pairs from GWAS1

ibd_gt_0.8_processed <- ibd_gt_0.8_processed %>%
  mutate(RemovedID = if_else(Remove == Batch1, IID1, IID2))

ibd_gt_0.8_processed[ibd_gt_0.8_processed$RemovedID %in% "A0405_Kamboh_C_0030", ]$RemovedID <- "A0406"
ibd_gt_0.8_processed[ibd_gt_0.8_processed$IID2 %in% "A0386_A0386", ]$Remove <- "GWAS1"
ibd_gt_0.8_processed[ibd_gt_0.8_processed$IID2 %in% "A0386_A0386", ]$RemovedID <- "A0386_A0386"

ibd_gt_0.8_rmList <- c(ibd_gt_0.8_processed[is.na(ibd_gt_0.8_processed$Keep2),]$IID1,
                       ibd_gt_0.8_processed[is.na(ibd_gt_0.8_processed$Keep2),]$IID2,
                       ibd_gt_0.8_processed[!is.na(ibd_gt_0.8_processed$Remove) & !is.na(ibd_gt_0.8_processed$Keep2),]$RemovedID, 
                       ibd_gt_0.8_processed[!is.na(ibd_gt_0.8_processed$Keep2) & is.na(ibd_gt_0.8_processed$Remove),]$IID1)
ibd_gt_0.8_rmList[duplicated(ibd_gt_0.8_rmList)] # 155
# "B0339_B0339"         "A1381_Kamboh_C_0086" "B0336_B0336"         "W0027_Kamboh_C_0594" "A0135_A0135"        
# "A1922_A1922"         "M0417_Kamboh_C_0347"
length(unique(ibd_gt_0.8_rmList))
# 148

### sheet 2
ibd_lt_0.8 <- read_excel("./QC_Check/IBD_Check_2024_06_10_Frank.xlsx", sheet = "PI_HAT < 0.8")

table(ibd_lt_0.8$Notes)
ibd_lt_0.8_rmList <- c(ibd_lt_0.8[!is.na(ibd_lt_0.8$Notes) & ibd_lt_0.8$Notes == "Identical IID and AlzID", ]$IID1, 
                       ibd_lt_0.8[!is.na(ibd_lt_0.8$Notes) & ibd_lt_0.8$Notes == "Identical IID and AlzID", ]$IID2)
ibd_lt_0.8_rmList[duplicated(ibd_lt_0.8_rmList)]
# a2058
length(unique(ibd_lt_0.8_rmList))
# 177

intersect(ibd_lt_0.8_rmList, ibd_gt_0.8_rmList)
# "A0208_A0208" "A2058_A2058" "S0113"      
### combine removal list
ibd_removal <- data.frame(V1 = c(ibd_lt_0.8_rmList, ibd_gt_0.8_rmList), 
                          V2 = c(ibd_lt_0.8_rmList, ibd_gt_0.8_rmList)) %>%
  distinct()
# 322

### assign Removal.Reasons
all_tmp <- all

for (seq in 1:nrow(all_tmp)){
  if(is.na(all_tmp[seq,"Removal.Reasons"])){
    all_tmp[seq, "Removal.Reasons"] <- ifelse(all_tmp[seq, "V1"] %in% ibd_removal$V1, "Cross-batch IBD Failed", NA)
  }
}

# add detailed ibd failure information
note1 <- "IBD PI_HAT >= 0.8 pair with conflicted Sex (in GWAS1, remove both)"
rm1 <- c(ibd_gt_0.8_processed[is.na(ibd_gt_0.8_processed$Keep2),]$IID1,
         ibd_gt_0.8_processed[is.na(ibd_gt_0.8_processed$Keep2),]$IID2)

note2 <- "IBD PI_HAT >= 0.8 pair with conflicted Sex (remove the one in GWAS1)" 
rm2 <- ibd_gt_0.8_processed[!is.na(ibd_gt_0.8_processed$Notes) & ibd_gt_0.8_processed$Notes == "Different IID and AlzID" & !is.na(ibd_gt_0.8_processed$Keep2) & ibd_gt_0.8_processed$Sex == 0,]$RemovedID

note3 <- "IBD PI_HAT >= 0.8 pair with different IID and AlzID (remove the one in the older batch)"
rm3 <- c(ibd_gt_0.8_processed[ibd_gt_0.8_processed$Sex == 1 & !is.na(ibd_gt_0.8_processed$Notes) & ibd_gt_0.8_processed$Notes == "Different IID and AlzID" & !is.na(ibd_gt_0.8_processed$RemovedID),]$RemovedID, 
         ibd_gt_0.8_processed[ibd_gt_0.8_processed$Sex == 1 & !is.na(ibd_gt_0.8_processed$Notes) & !is.na(ibd_gt_0.8_processed$Keep2) & is.na(ibd_gt_0.8_processed$Remove) & ibd_gt_0.8_processed$Notes == "Different IID and AlzID", ]$IID1)

note4 <-  "IBD PI_HAT >= 0.8 pair with different IID but identical AlzID (remove the one in the older batch)"
rm4 <- c(ibd_gt_0.8_processed[ibd_gt_0.8_processed$Sex == 1 & !is.na(ibd_gt_0.8_processed$Notes) & ibd_gt_0.8_processed$Notes == "Different IID but Identical AlzID" & !is.na(ibd_gt_0.8_processed$RemovedID),]$RemovedID, 
         ibd_gt_0.8_processed[ibd_gt_0.8_processed$Sex == 1 & !is.na(ibd_gt_0.8_processed$Notes) & !is.na(ibd_gt_0.8_processed$Keep2) & is.na(ibd_gt_0.8_processed$Remove) & ibd_gt_0.8_processed$Notes == "Different IID but Identical AlzID", ]$IID1)

note5 <- "IBD PI_HAT >= 0.8 pair with identical IID and AlzID (remove the one in the older batch)"
rm5 <- c(ibd_gt_0.8_processed[ibd_gt_0.8_processed$Sex == 1 & !is.na(ibd_gt_0.8_processed$Notes) & ibd_gt_0.8_processed$Notes == "Identical IID and AlzID" & !is.na(ibd_gt_0.8_processed$RemovedID),]$RemovedID, 
         ibd_gt_0.8_processed[ibd_gt_0.8_processed$Sex == 1 & !is.na(ibd_gt_0.8_processed$Notes) & !is.na(ibd_gt_0.8_processed$Keep2) & is.na(ibd_gt_0.8_processed$Remove) & ibd_gt_0.8_processed$Notes == "Identical IID and AlzID", ]$IID1)

note6 <- "IBD PI_HAT < 0.8 pair with identical IID and AlzID (remove both)"
rm6 <- c(ibd_lt_0.8[!is.na(ibd_lt_0.8$Notes) & ibd_lt_0.8$Notes == "Identical IID and AlzID", ]$IID1, 
         ibd_lt_0.8[!is.na(ibd_lt_0.8$Notes) & ibd_lt_0.8$Notes == "Identical IID and AlzID", ]$IID2)

# rm_check <- c(setdiff(c(rm1,rm2,rm3,rm4,rm5,rm6), ibd_removal$V1),
#               setdiff(ibd_removal$V1, c(rm1,rm2,rm3,rm4,rm5,rm6)))
# empty 

# combine notes and rm lists into a named list
rm_notes <- list(note1 = rm1, note2 = rm2, note3 = rm3, note4 = rm4, note5 = rm5, note6 = rm6)

# function to accumulate notes based on presence in rm lists
get_notes <- function(V1, rm_notes) {
  notes <- sapply(names(rm_notes), function(n) ifelse(V1 %in% rm_notes[[n]], n, NA))
  paste(na.omit(notes), collapse = ";")
}

# Add Notes column based on the conditions
ibd_removal_tmp <- ibd_removal %>%
  rowwise() %>%
  mutate(Notes = get_notes(V1, rm_notes)) %>%
  mutate(Notes.full = get_notes(V1, rm_notes)) %>%
  ungroup()

ibd_removal_tmp$Notes.full <- gsub("note1", note1, ibd_removal_tmp$Notes.full)
ibd_removal_tmp$Notes.full <- gsub("note2", note2, ibd_removal_tmp$Notes.full)
ibd_removal_tmp$Notes.full <- gsub("note3", note3, ibd_removal_tmp$Notes.full)
ibd_removal_tmp$Notes.full <- gsub("note4", note4, ibd_removal_tmp$Notes.full)
ibd_removal_tmp$Notes.full <- gsub("note5", note5, ibd_removal_tmp$Notes.full)
ibd_removal_tmp$Notes.full <- gsub("note6", note6, ibd_removal_tmp$Notes.full)

write.table(ibd_removal_tmp[c("V1","V2")], "./QC_Check/IBD_Remove_list.txt", col.names = FALSE, row.names = FALSE, sep = " ", quote = FALSE)
write.table(ibd_removal_tmp, "./QC_Check/IBD_Remove_list_with_detailed_reasons.txt", col.names = FALSE, row.names = FALSE, sep = " ", quote = FALSE)

###--------------------------------------------------------------------------###
### IBD double check
###--------------------------------------------------------------------------###

ibd_double_check <- fread("./QC_Check/AllCohorts_AllChr_IBD2.genome")
dim(ibd_double_check)
# 323  14
nrow(ibd_double_check[ibd_double_check$PI_HAT >= 0.8,])

ibd_double_check1 <- merge(ibd_double_check[,c("FID1", "IID1", "FID2", "IID2", "PI_HAT")], 
                           all[all$GWAS_QC == "Pass",c("V1","AlzID","Study","Batch","Sex","Race")], 
                           by.x = "IID1", by.y = "V1", all.x = TRUE) %>% distinct()

colnames(ibd_double_check1) <- c("IID1","FID1","FID2","IID2","PI_HAT","AlzID1","Study1","Batch1","Sex1","Race1")

ibd_double_check2 <- merge(ibd_double_check1, 
                           all[all$GWAS_QC == "Pass",c("V1","AlzID","Study","Batch","Sex","Race")], 
                           by.x = "IID2", by.y = "V1", all.x = TRUE) %>% distinct()

colnames(ibd_double_check2) <- c("IID2","IID1","FID1","FID2","PI_HAT","AlzID1","Study1","Batch1","Sex1","Race1",
                                 "AlzID2","Study2","Batch2","Sex2","Race2")

ibd_check2 <- ibd_double_check2[order(ibd_double_check2$PI_HAT, decreasing = TRUE), 
                                c("PI_HAT","FID1","IID1","AlzID1","Study1","Batch1","Sex1","Race1",
                                  "FID2","IID2","AlzID2","Study2","Batch2","Sex2","Race2")]

ibd_check2$Notes <- as.character(NA_character_)
for (seq in 1:nrow(ibd_check2)){
  ibd_check2[seq, "Notes"] <- ifelse(grepl(toupper(ibd_check2[seq,"IID1"]), toupper(ibd_check2[seq,"IID2"])) | grepl(toupper(ibd_check2[seq,"IID2"]), toupper(ibd_check2[seq,"IID1"])) & 
                                       ibd_check2[seq,"AlzID1"] %in% ibd_check2[seq,"AlzID2"], "Identical IID and AlzID", 
                                     ifelse(ibd_check2[seq,"AlzID1"] %in%  ibd_check2[seq,"AlzID2"], "Different IID but Identical AlzID", 
                                            ifelse(!(ibd_check2[seq,"AlzID1"] %in% ibd_check2[seq,"AlzID2"]) & ibd_check2[seq,"PI_HAT"] > 0.8, "Different IID and AlzID", NA)))
}

setdiff(ibd_check[ibd_check$Notes == "Different IID but Identical AlzID" & ibd_check$PI_HAT < 0.8, ]$AlzID1, ibd_check2[ibd_check2$Notes == "Different IID but Identical AlzID", ]$AlzID1)
# "ALZ100818" "ALZ105507"

###--------------------------------------------------------------------------###
### prepare two lists that may need modification 
###--------------------------------------------------------------------------###

# list 1: Pairs that genotype indicates identical), but different AlzIDs
list1 <- ibd_gt_0.8_processed[!is.na(ibd_gt_0.8_processed$Notes) & ibd_gt_0.8_processed$Notes == "Different IID and AlzID", 
                              c("PI_HAT","IID1","AlzID1","Study1","IID2","AlzID2","Study2")]

# list 2: Pairs with identical AlzIDs, but different LabIDs (genotype indicates differences)
list2 <- ibd_lt_0.8[!is.na(ibd_lt_0.8$Notes) & ibd_lt_0.8$Notes == "Different IID but Identical AlzID", 
                              c("PI_HAT","IID1","AlzID1","Study1","IID2","AlzID2","Study2")]

wb <- createWorkbook()

addWorksheet(wb, "Sheet1")
writeData(wb, "Sheet1", list1)

addWorksheet(wb, "Sheet2")
writeData(wb, "Sheet2", list2)

saveWorkbook(wb, "./QC_Check/after_IBD_Check_ID_Discrepancies.xlsx", overwrite = TRUE)

IBD_failed <- all_tmp[all_tmp$Removal.Reasons %in% "Cross-batch IBD Failed",]
write.xlsx(IBD_failed[order(IBD_failed$Batch),], "./QC_Check/GWAS_Harmonization_IBD_Failed.xlsx")

###--------------------------------------------------------------------------###
### update the final list
###--------------------------------------------------------------------------###

# update the final list
all_tmp$Additional_Notes <- as.character(NA_character_)
for (seq in 1:nrow(all_tmp)){
  scanName <- all_tmp[seq, "V1"]
  if(!is.na(all_tmp[seq, "Removal.Reasons"]) & all_tmp[seq, "Removal.Reasons"] == "Cross-batch IBD Failed"){
    all_tmp[seq,"Additional_Notes"] <- ifelse(scanName %in% ibd_removal_tmp$V1, ibd_removal_tmp[ibd_removal_tmp$V1 == scanName, "Notes.full"], NA)
  }
}

table(all_tmp$Batch, all_tmp$Removal.Reasons)

# update all_tmp GWAS_QC
all_tmp$GWAS_QC <- ifelse(is.na(all_tmp$Removal.Reasons), "Pass", "Failed")
table(all_tmp$GWAS_QC)
# Failed   Pass 
#    956  10431

table(all_tmp[all_tmp$GWAS_QC == "Pass",]$Batch)

# manually assign race to BENIN and SLVDS group
all_tmp[all_tmp$Stud == "BENIN", "Race"] <- "Black"
all_tmp[all_tmp$Stud == "SLVDS", "Race"] <- "White"

# for those marked as Different IID but Identical AlzID, assign -12 to these samples
all_tmp[all_tmp$V1 %in% c(list2$IID1, list2$IID2), "AlzID"] <- "-12"
# length(unique(c(list2$IID1, list2$IID2)))

write.xlsx(all_tmp, "./QC_Check/GWAS_Harmonization_N=11387_Demographic_Final.xlsx")

