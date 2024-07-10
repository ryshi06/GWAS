
rm(list = ls())
gc()

library(ggplot2)
library(data.table)
library(dplyr)
library(openxlsx)

setwd("/ix/kfan/Ruyu/harmonization/All_Cohorts")

# rs429358
# 19:45411941:T:C[b37]_C

# rs141622900 (in LD with rs7412) 
# 19:45426792:G:A[b37]_A
# correlated alleles: C=G,T=A	

raw <- fread("./QC_Check/AllCohorts_AllChr_Annot_APOE_SNPs.raw")
ped <- fread("./QC_Check/AllCohorts_AllChr_Annot_APOE_SNPs.ped", sep = " ")

raw$apoe <- ifelse(raw$`19:45411941:T:C[b37]_C` == 0 & raw$`19:45426792:G:A[b37]_A` == 2, "22",
                   ifelse(raw$`19:45411941:T:C[b37]_C` == 0 & raw$`19:45426792:G:A[b37]_A` == 1, "23", 
                          ifelse(raw$`19:45411941:T:C[b37]_C` == 1 & raw$`19:45426792:G:A[b37]_A` == 1, "24",
                                 ifelse(raw$`19:45411941:T:C[b37]_C` == 0 & raw$`19:45426792:G:A[b37]_A` == 0, "33",
                                        ifelse(raw$`19:45411941:T:C[b37]_C` == 1 & raw$`19:45426792:G:A[b37]_A` == 0, "34", "44")))))

table(raw$apoe)
# 22   23   24   33   34   44 
# 43  773  210 6038 2919  448 

###--------------------------------------------------------------------------###
### load apoe phenotype
###--------------------------------------------------------------------------###

pheno <- read.xlsx("./Harmonization_Phenotype_2024_06_15.xlsx")

length(intersect(raw$IID, pheno$scanName))
dd <- merge(raw[,c("IID","19:45411941:T:C[b37]_C","19:45426792:G:A[b37]_A","apoe")], pheno[pheno$GWAS_QC == "Pass", ], by.x = "IID", by.y = "scanName") %>% distinct()

table(dd$APOE)
# 22   23   24   33   34   44 
# 49 1043  247 5605 2796  452 
sum(is.na(dd$APOE))
# 239

###--------------------------------------------------------------------------###
### check APOE mismatch
###--------------------------------------------------------------------------###

### subset to unmatched individuals
mismatch <- dd[!is.na(dd$apoe) & !is.na(dd$APOE) & dd$apoe != dd$APOE,]
dim(mismatch)
# 1002   17
mismatch$alert <- as.character(paste0(as.character(mismatch$apoe), "|", as.character(mismatch$APOE)))
mismatch$alert <- factor(mismatch$alert, levels = c("22|23","22|33",
                                                    "23|22","23|24","23|33","23|34",
                                                    "24|23","24|33","24|34",
                                                    "33|22","33|23","33|24","33|34","33|44",
                                                    "34|22","34|23","34|24","34|33","34|44",
                                                    "44|23","44|33","44|34"))
table(mismatch$alert)
table(mismatch$alert, mismatch$Batch)

no_apoe <- dd[is.na(dd$APOE),]
dim(no_apoe)
# 265  18
table(no_apoe$Batch)
table(no_apoe$Batch, no_apoe$StudyList)

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
  
}

table(dd$apoe_alert)
# 22|23               22|33               23|22               23|24 
#     7                  19                  22                   1 
# 23|33               23|34               24|23               24|33 
#    81                  10                   8                  26 
# 24|34               33|22               33|23               33|24 
#    23                   8                 362                   9 
# 33|34               33|44               34|22               34|23 
#   133                  16                   2                  24 
# 34|24               34|33               34|44               44|23 
#    87                 106                  22                   1 
# 44|33               44|34 no APOE in database 
#    20                  10                 265

table(dd$Batch, dd$apoe_alert)
d1 <- dd[!is.na(dd$apoe_alert) & dd$apoe_alert == "no APOE in database",]
write.xlsx(d1[order(d1$Batch),], "./QC_Check/AllCohorts_AllChr_Harmonized_no_APOE.xlsx")

d2 <- dd[!is.na(dd$apoe_alert) & dd$apoe_alert != "no APOE in database",]
write.xlsx(d2[order(d2$Batch),], "./QC_Check/AllCohorts_AllChr_Harmonized_APOE_discrepancies.xlsx")

###--------------------------------------------------------------------------###
### update APOE genotype based on database
###--------------------------------------------------------------------------###

head(ped)
ped2 <- dd[,c("IID","apoe", "APOE", "apoe_alert")]

ped2$apoe_genotype <- ifelse(dd$apoe_alert %in% "no APOE in database", dd$apoe, dd$APOE)
table(ped2$apoe_genotype)
# 22   23   24   33   34   44 
# 49 1064  248 5764 2852  454
sum(is.na(ped2$apoe_genotype))

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

##write.table(ped3, "./QC_Check/AllCohorts_AllChr_Annot_APOE_SNPs.ped", sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)

# map <- fread("./QC_Check/AllCohorts_AllChr_Annot_APOE_SNPs.map")

###--------------------------------------------------------------------------###
### after-merge APOE SNP check
###--------------------------------------------------------------------------###

new.raw <- fread("/ix/kfan/Ruyu/harmonization/All_Cohorts/QC_Check/AllCohorts_AllChr_Annot_APOE_SNPs_updated.raw")

new.raw$apoe <- ifelse(new.raw$`19:45411941:T:C[b37]_C` == 0 & new.raw$`19:45426792:G:A[b37]_A` == 2, "22",
                       ifelse(new.raw$`19:45411941:T:C[b37]_C` == 0 & new.raw$`19:45426792:G:A[b37]_A` == 1, "23", 
                              ifelse(new.raw$`19:45411941:T:C[b37]_C` == 1 & new.raw$`19:45426792:G:A[b37]_A` == 1, "24",
                                     ifelse(new.raw$`19:45411941:T:C[b37]_C` == 0 & new.raw$`19:45426792:G:A[b37]_A` == 0, "33",
                                            ifelse(new.raw$`19:45411941:T:C[b37]_C` == 1 & new.raw$`19:45426792:G:A[b37]_A` == 0, "34", "44")))))

table(ped2$apoe_genotype)
# 22   23   24   33   34   44 
# 49 1061  250 5762 2854  455

table(new.raw$apoe)
# 22   23   24   33   34   44 
# 49 1061  250 5762 2854  455

# finished

