
rm(list = ls())
gc()

library(ggplot2)
library(data.table)
library(dplyr)

setwd("/ix/kfan/Ruyu/harmonization/All_Cohorts")

# rs429358
# 19:45411941:T:C[b37]_C

# rs141622900 (in LD with rs7412) 
# 19:45426792:G:A[b37]_A
# correlated alleles: C=G,T=A	

raw <- fread("./AllCohorts_AllChr_Annot_APOE_SNPs.raw")
ped <- fread("./AllCohorts_AllChr_Annot_APOE_SNPs.ped")

raw$apoe <- ifelse(raw$`19:45411941:T:C[b37]_C` == 0 & raw$`19:45426792:G:A[b37]_A` == 2, "22",
                   ifelse(raw$`19:45411941:T:C[b37]_C` == 0 & raw$`19:45426792:G:A[b37]_A` == 1, "23", 
                          ifelse(raw$`19:45411941:T:C[b37]_C` == 1 & raw$`19:45426792:G:A[b37]_A` == 1, "24",
                                 ifelse(raw$`19:45411941:T:C[b37]_C` == 0 & raw$`19:45426792:G:A[b37]_A` == 0, "33",
                                        ifelse(raw$`19:45411941:T:C[b37]_C` == 1 & raw$`19:45426792:G:A[b37]_A` == 0, "34", "44")))))

table(raw$apoe)
# 22   23   24   33   34   44 
# 44  779  219 6107 2999  463 

###--------------------------------------------------------------------------###
### load apoe phenotype
###--------------------------------------------------------------------------###

d <- fread("./All_Cohorts_2nd_Imputation_Phenotype.txt")

length(intersect(raw$IID, d$scanName))
dd <- merge(raw[,c("IID","19:45411941:T:C[b37]_C","19:45426792:G:A[b37]_A","apoe")], d, by.x = "IID", by.y = "scanName")

table(dd$apoe_recode)
# 22   23   24   33   34   44 
# 49 1045  252 5645 2841  457
sum(is.na(dd$apoe_recode))
# 313

###--------------------------------------------------------------------------###
### check APOE mismatch
###--------------------------------------------------------------------------###

### subset to unmatched individuals
mismatch <- dd[!is.na(dd$apoe) & !is.na(dd$apoe_recode) & dd$apoe != dd$apoe_recode,]
dim(mismatch)
# 1038   10
mismatch$alert <- as.character(paste0(as.character(mismatch$apoe), "|", as.character(mismatch$apoe_recode)))
mismatch$alert <- factor(mismatch$alert, levels = c("22|23","22|33",
                                                    "23|22","23|24","23|33","23|34",
                                                    "24|23","24|33","24|34",
                                                    "33|22","33|23","33|24","33|34","33|44",
                                                    "34|22","34|23","34|24","34|33","34|44",
                                                    "44|23","44|33","44|34"))
table(mismatch$alert)
table(mismatch$alert, mismatch$Batch)

no_apoe <- dd[is.na(dd$apoe_recode),]
dim(no_apoe)
# 313  10

# race discrepancies
race_mismatch <- fread("./QC_Check/AllCohorts_AllChr_Race_Discrepancies.txt")
table(race_mismatch$alert)
table(race_mismatch$alert, race_mismatch$Batch)

no_race <- dd[dd$race_recode == "Unknown",]

###--------------------------------------------------------------------------###
### assign back apoe/race alert
###--------------------------------------------------------------------------###

for (seq in 1:nrow(dd)){
  
  id <- dd[seq,]$IID
  
  # check apoe discrepancies
  tmp1 <- mismatch[mismatch$IID == id, ]
  tmp2 <- no_apoe[no_apoe$IID == id,]
  if(nrow(tmp1) > 0){
    dd$apoe_alert[seq] <- as.character(tmp1$alert)
  } else if (nrow(tmp2) > 0){
    dd$apoe_alert[seq] <- "no APOE in database"
  } else {
    dd$apoe_alert[seq] <- NA
  }
  
  # check race discrepancies
  tmp1 <- race_mismatch[race_mismatch$IID == id, ]
  tmp2 <- no_race[no_race$IID == id , ]
  if(nrow(tmp1) > 0){
    dd$race_alert[seq] <- tmp1$alert
  } else if (nrow(tmp2) > 0){
    dd$race_alert[seq] <- "no Race in database"
  } else {
    dd$race_alert[seq] <- NA
  }
  
}

# check if the assignment is correct
table(dd$race_alert)
# no Race in database self-reported Asian clustered in White 
#                1340                                      6 
# self-reported Black clustered in White self-reported White clustered in Black 
#                                      8                                     25 
sum(is.na(dd$race_alert))
# 9232 (1379)
table(dd$race_alert, dd$Batch)

table(dd$apoe_alert)
# 22|23               22|33               23|22               23|24               23|33               23|34 
#     7                  20                  21                   1                  86                  10 
# 24|23               24|33               24|34               33|22               33|23               33|24 
#    10                  28                  24                   8                 364                   9 
# 33|34               33|44               34|22               34|23               34|24               34|33 
#   141                  21                   3                  25                  88                 113 
# 34|44               44|23               44|33               44|34 no APOE in database 
#    22                   1                  26                  10                 313

##openxlsx::write.xlsx(dd, "./QC_Check/AllCohorts_AllChr_Harmonized_Race_APOE_discrepancies.xlsx")

###--------------------------------------------------------------------------###
### update APOE genotype based on database
###--------------------------------------------------------------------------###

head(ped)
ped2 <- dd[,c("IID","apoe", "apoe_recode", "apoe_alert")]

ped2$apoe_genotype <- ifelse(!is.na(dd$apoe_alert) & dd$apoe_alert == "no APOE in database", dd$apoe, dd$apoe_recode)

ped2$`19:45411941:T:C[b37]` <- ifelse(ped2$apoe_genotype == "22", "TT", 
                                      ifelse(ped2$apoe_genotype == "23", "TT", 
                                             ifelse(ped2$apoe_genotype == "24", "CT", 
                                                    ifelse(ped2$apoe_genotype == "33", "TT", 
                                                           ifelse(ped2$apoe_genotype == "34", "CT", "CC")))))
ped2$`19:45426792:G:A[b37]` <- ifelse(ped2$apoe_genotype == "22", "AA", 
                                      ifelse(ped2$apoe_genotype == "23", "AG", 
                                             ifelse(ped2$apoe_genotype == "24", "AG", 
                                                    ifelse(ped2$apoe_genotype == "33", "GG", 
                                                           ifelse(ped2$apoe_genotype == "34", "GG", "GG")))))

split_digits <- function(x) {
  return(unlist(strsplit(as.character(x), "")))
}

for (seq in 1:nrow(ped2)){
  
  ped2$`19:45411941:T:C[b37].1`[seq] <- split_digits(ped2$`19:45411941:T:C[b37]`[seq])[1]
  ped2$`19:45411941:T:C[b37].2`[seq] <- split_digits(ped2$`19:45411941:T:C[b37]`[seq])[2]
  
  ped2$`19:45426792:G:A[b37].1`[seq] <- split_digits(ped2$`19:45426792:G:A[b37]`[seq])[1]
  ped2$`19:45426792:G:A[b37].2`[seq] <- split_digits(ped2$`19:45426792:G:A[b37]`[seq])[2]
  
}

ped_tmp <- merge(ped, ped2[,c("IID", "19:45411941:T:C[b37].1", "19:45411941:T:C[b37].2", "19:45426792:G:A[b37].1", "19:45426792:G:A[b37].2", "apoe_alert")], by.x = "V1", by.y = "IID")

ped_tmp$alert1 <- ifelse(ped_tmp$V7 != ped_tmp$`19:45411941:T:C[b37].1`, "Y", NA)
ped_tmp$alert2 <- ifelse(ped_tmp$V8 != ped_tmp$`19:45411941:T:C[b37].2`, "Y", NA)
ped_tmp$alert3 <- ifelse(ped_tmp$V9 != ped_tmp$`19:45426792:G:A[b37].1`, "Y", NA)
ped_tmp$alert4 <- ifelse(ped_tmp$V10 != ped_tmp$`19:45426792:G:A[b37].2`, "Y", NA)

tmp <- ped_tmp[!is.na(ped_tmp$alert1) | !is.na(ped_tmp$alert2) | !is.na(ped_tmp$alert2) | !is.na(ped_tmp$alert2) | !is.na(ped_tmp$apoe_alert) & ped_tmp$apoe_alert != "no APOE in database", ]

ped3 <- merge(ped[,paste0("V",1:6)], ped2[,c("IID","19:45411941:T:C[b37].1","19:45411941:T:C[b37].2","19:45426792:G:A[b37].1","19:45426792:G:A[b37].2")],
              by.x = "V1", by.y = "IID")
colnames(ped3) <- paste0("V",1:10)

ped3 <- ped3[match(ped$V1, ped3$V1), ]
all.equal(ped3$V1, ped$V1)

write.table(ped3, "./AllCohorts_AllChr_Annot_APOE_SNPs.ped", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# map <- fread("./QC_Check/AllCohorts_AllChr_Annot_APOE_SNPs.map")

###--------------------------------------------------------------------------###
### after-merge APOE SNP check
###--------------------------------------------------------------------------###

new.raw <- fread("/ix/kfan/Ruyu/harmonization/All_Cohorts/AllCohorts_AllChr_Annot_APOE_SNPs_updated.raw")

new.raw$apoe <- ifelse(new.raw$`19:45411941:T:C[b37]_C` == 0 & new.raw$`19:45426792:G:A[b37]_A` == 2, "22",
                       ifelse(new.raw$`19:45411941:T:C[b37]_C` == 0 & new.raw$`19:45426792:G:A[b37]_A` == 1, "23", 
                              ifelse(new.raw$`19:45411941:T:C[b37]_C` == 1 & new.raw$`19:45426792:G:A[b37]_A` == 1, "24",
                                     ifelse(new.raw$`19:45411941:T:C[b37]_C` == 0 & new.raw$`19:45426792:G:A[b37]_A` == 0, "33",
                                            ifelse(new.raw$`19:45411941:T:C[b37]_C` == 1 & new.raw$`19:45426792:G:A[b37]_A` == 0, "34", "44")))))

table(ped2$apoe_genotype)
# 22   23   24   33   34   44 
# 49 1068  255 5837 2933  469 

table(new.raw$apoe)
# 22   23   24   33   34   44 
# 49 1067  255 5838 2933  469 

# finished

