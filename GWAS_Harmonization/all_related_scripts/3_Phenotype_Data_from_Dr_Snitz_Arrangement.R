
rm(list=ls())
setwd("/home/Data/GEM_GWAS/Frank")
options(stringsAsFactors = FALSE)
library(foreign)
library(gdata)
library(haven)

###--------------------------------------------------------------------------###
###
### This script is trying to collect phenoytpes that 
###   we are going to use in the GEM GWAS study
###
###                                       Frank
###                                       07/09/2109
###
###--------------------------------------------------------------------------###

###--------------------------------------------------------------------------###
### Notes from Dr. Snitz's email
###--------------------------------------------------------------------------###

# (1) GEMS endpoint status (status at end of the study):
#         identifies incident dementia cases. 
#         Also has date of birth so age can be calculated at any follow-up date.
# (2) Serial neuropsychological test data: 
#         for domain score calculation. 
#         I will have to send you a list (soon) of tests for constructing 
#         domain composite scores – we can try to harmonize with MYHAT as much as possible.
# (3) APOE
# (4) Data dictionary: Excel file 
# (5) Exit interview:
#         which has date of close-out visit, 
#         can be used to calculate age at study close-out in 2007-2008.
# (6) Serial CDR scores:
#         (variable of interest = ‘cdglobl’, which is the global CDR rating)
# (7) Serial ADAS-cog: 
#         cognitive screening test done every 6 months
# (8) Serial 3MSE: 
#         cognitive screening test done every 6 months

###--------------------------------------------------------------------------###
### BASIC INFORMATION (Not Logitudinal Type Data)
###--------------------------------------------------------------------------###

# basic info

base <- read.spss("../Public/Useful_Data/GEMBAS~1.sav", to.data.frame = TRUE)

# ID, sex, race, age of enrollment, baseline cdr, clinic,
# education level, height, weight, BMI, 

d <- subset(base, select=c("idno","gender","race","agernd", "edu", "clinic",
                           "cdr1", "vsheight", "vsweight", "bmi"))

d$GEM_ID <- d$idno

d$sex <- ifelse(d$gender=="0: Female", "F", "M")
# table(d$gender, d$sex, useNA = "ifany")  

d$race2 <- gsub(" ","", as.character(d$race))
d$race <- sapply(d$race2, function(x){
                  unlist(strsplit(x, ":"))[2]})
# table(d$race2, d$race, useNA = "ifany")

d$age <- d$agernd

d$center <- sapply(as.character(d$clinic), function(x){
                    unlist(strsplit(x, ": "))[2]})
# table(d$center, d$clinic, useNA="ifany")

d$cdr <- d$cdr1
d$weight <- d$vsweight
d$height <- d$vsheight

GEM_BASE <- subset(d, select=c("GEM_ID", "sex", "race", "age", "edu", "center",
                               "cdr", "weight", "height", "bmi"))

str(GEM_BASE)
rm(list=setdiff(ls(),"GEM_BASE"))

###--------------------------------------------------------------------------###

# Other Basic Information not in the baseline

# exit_interview <- read_sav("exitinterview.sav")
# str(exit_interview)
# only has 2025 records

end_points <- read_sav("../Public/Useful_Data/GEMendpoints.sav")
# str(end_points)
# missing 1 sample

# str(end_points$dementia)
# table(end_points$dementia, useNA="ifany")

d2 <- subset(end_points, select=c("idno","dementia"))
d2$dementia[is.na(d2$dementia)] <- 3

d2$AD <- ifelse(d2$dementia==1, "AD", "Control")
d2$Dementia <- ifelse(d2$dementia==1, "AD", 
               ifelse(d2$dementia==2, "Non-AD", "Control"))

# str(d2)
# table(d2$dementia, d2$AD)
# table(d2$dementia, d2$Dementia)

names(d2) <- c("GEM_ID", "dementia_old", "AD","dementia")
d3 <- as.data.frame(subset(d2, select = c("GEM_ID","AD","dementia")))

GEM_CLIN <- merge(GEM_BASE, d3, all=TRUE)


apoe <- read.xls("../Public/Useful_Data/APOE_GEMS.xls", sheet = 1)
str(apoe)
# 2455 (maybe this is the data that we gave to Dr. Snitz)
# We need to merge this file later with everyone's APOE genotype

save(GEM_CLIN, file="Basic_Clinical_Info_for_GEM.RData")


# str(GEM_CLIN)

###--------------------------------------------------------------------------###
### Logitudinal Type Data (No Run !!)       8/22/2019 Frank
###--------------------------------------------------------------------------###

# # Serial neuropsychological test data
# neuro <- read.spss("neuropsych.sav", to.data.frame = TRUE)
# str(neuro)
# 
# # MMSE, CDR, ADAS
# X3MSE <- read_sav("3MSE.sav")
# str(X3MSE)
# adas <- read.spss("adas.sav", to.data.frame = TRUE)
# str(adas)
# cdr <- read.spss("cdr.sav", to.data.frame = TRUE)
# str(cdr)





