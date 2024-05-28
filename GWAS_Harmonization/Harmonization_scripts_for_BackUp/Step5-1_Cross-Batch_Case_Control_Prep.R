
setwd("/ix/kfan/Ruyu/harmonization/All_Cohorts/")

pattern <- "_allChr_scanName.txt"

files <- list.files(path = "/ix/kfan/Ruyu/harmonization/All_Cohorts/pseudoGWAS_input", pattern = "_allChr_scanName.txt")

scanName <- list()
for (file in files){
  scanName[[file]] <- read.table(paste0("/ix/kfan/Ruyu/harmonization/All_Cohorts/pseudoGWAS_input/",file))
  scanName[[file]]$scanName <- stringr::str_sub(scanName[[file]]$V1, 1, nchar(scanName[[file]]$V1)/2)
  scanName[[file]]$cohort <- gsub("_allChr_scanName.txt","",file)
}

sample <- do.call(rbind, scanName)
table(sample$cohort)
dim(sample)

### -------------------------------------------------------------------------###
### load data from Lynda
### -------------------------------------------------------------------------###

dataRequest <- openxlsx::read.xlsx("V_DataRequest_Ruyu_Harmonization_2024_05_01.xlsx")
dataRequest$ori_Batches <- NA
dim(dataRequest)
for (seq in 1:nrow(dataRequest)){
  scan_name <- dataRequest[seq,]$Scan_name
  batch <- unique(sample[sample$scanName == scan_name,]$cohort)
  if(length(batch)>1){
    dataRequest[dataRequest$Scan_name == scan_name,]$ori_Batches <- paste(batch, collapse = ",")
  } else {
    dataRequest[dataRequest$Scan_name == scan_name,]$ori_Batches <- batch
  }
}

# View(dataRequest[is.na(dataRequest$Case_Control_Current),])
# View(dataRequest[is.na(dataRequest$GWAS_Batches),])
# nrow(dataRequest[is.na(dataRequest$Case_Control_Current),])
# table(dataRequest[is.na(dataRequest$Case_Control_Current),]$ori_Batches)
# table(dataRequest[is.na(dataRequest$Case_Control_Current),]$Study)

# assign those with missing Case_Control_Current to Cases
dataRequest[is.na(dataRequest$Case_Control_Current),]$Case_Control_Current <- "Case"
dataRequest[dataRequest$Case_Control_Current == -1 & dataRequest$Study == "ADRC",]$Case_Control_Current <- "Case"
dataRequest[dataRequest$Case_Control_Current == -1 & dataRequest$Study %in% c("GEM","MOVIES"),]$Case_Control_Current <- "Control"
dataRequest[dataRequest$Case_Control_Current == "EOAD",]$Case_Control_Current <- "Case"

# for those duplicated in cohorts, assign them to the cohort with less samples
table(sample$cohort)

gwas1_size <- 2395
gwas2_size <- 711
gwas3_size <- 613
gwas4_size <- 4318
gwas5_size <- 245
gem_size <- 2738

dataRequest$final_Batch <- NA
for (seq in 1:nrow(dataRequest)){
  
  scan_name <- dataRequest[seq,]$Scan_name
  old_batch <- dataRequest[dataRequest$Scan_name == scan_name, ]$ori_Batches
  
  if(length(grep(",", old_batch)) > 0){
    batches <- strsplit(old_batch, ",")[[1]]
    new_batch <- ifelse(length(grep("GWAS5", batches)) > 0, "GWAS5", 
                        ifelse(length(grep("GWAS3", batches)) > 0, "GWAS3", 
                               ifelse(length(grep("GWAS2", batches)) > 0, "GWAS2", 
                                      ifelse(length(grep("GWAS1", batches)) > 0, "GWAS1", 
                                             ifelse(length(grep("GEM", batches)) > 0, "GEM", "GWAS4")))))
  } else {
    new_batch <- old_batch
  }
  
  dataRequest[dataRequest$Scan_name == scan_name, ]$final_Batch <- new_batch
}

table(dataRequest$final_Batch)
sum(is.na(dataRequest$final_Batch))
table(dataRequest$Case_Control_Current)
sum(is.na(dataRequest$Case_Control_Current))

# for those controls without Sex information, they are from GWAS4, and their biological sexes are Female
dataRequest[dataRequest$Sex %in% c(-1,-11),]$Sex <- "F"
# sum(is.na(dataRequest$Sex))

### merge the fam scanName list with the dataRequest
dataRequest$V1 <- paste0(dataRequest$Scan_name,"_",dataRequest$Scan_name)

fam <- fread("./pseudoGWAS_input/All_Cohorts_allChr_noDups.fam")
sample <- as.data.frame(unique(sub(".*:", "", fam$V1)))
colnames(sample) <- "V1"
dd <- merge(sample, dataRequest[,c("V1", "Sex","Race","Case_Control_Current", "final_Batch")], by = "V1", all.x = TRUE)

table(dd$final_Batch)

# all samples not in dataRequest are from GWAS4, all Controls, all Males (those samples were excluded before due to manually assigned "F" and marked as mislabeled sex)

dd[is.na(dd$final_Batch),]$Case_Control_Current <- "Control"
dd[is.na(dd$final_Batch),]$Sex <- "M"
dd[is.na(dd$final_Batch),]$final_Batch <- "GWAS4"

dd_ids <- as.data.frame(cbind(dd$V1, dd$V1))
# write.table(dd_ids, "./pseudoGWAS_input/plink_IDs_noDup.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
### last check
dd <- dd %>% distinct()
table(dd$Case_Control_Current)
# Case Control 
# 3261    7555 
sum(is.na(dd$Case_Control_Current))
# 0 

### generate case/control list for each cohort

cases <- as.data.frame(cbind(dd[dd$Case_Control_Current == "Case", ]$V1,
                             dd[dd$Case_Control_Current == "Case", ]$V1))
colnames(cases) <- c("FID", "IID")
cases <- cases %>% distinct()
dim(cases)
# [1] 3261    2
write.table(cases, "./pseudoGWAS_input/AllCohorts_Cases_ID_forPlink.txt", 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

controls <- as.data.frame(cbind(dd[dd$Case_Control_Current == "Control", ]$V1,
                                dd[dd$Case_Control_Current == "Control", ]$V1))
colnames(controls) <- c("FID", "IID")
controls <- controls %>% distinct()
dim(controls)
# [1] 7555    2
# write.table(controls, "./pseudoGWAS_input/AllCohorts_Controls_ID_forPlink.txt", 
#             col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

### prepare phenotype file for pseudoGWAS

control_d <- dd[dd$Case_Control_Current == "Control", c("V1", "final_Batch", "Sex", "Race")]
table(control_d$final_Batch)

# PLINK assumes that controls are coded as 1 and cases are coded as 2 by default
control_d$GWAS1_2 <- ifelse(control_d$final_Batch == "GWAS1", 2, 
                            ifelse(control_d$final_Batch == "GWAS2", 1, NA))
control_d$GWAS1_3 <- ifelse(control_d$final_Batch == "GWAS1", 2, 
                            ifelse(control_d$final_Batch == "GWAS3", 1, NA))
control_d$GWAS1_4 <- ifelse(control_d$final_Batch == "GWAS1", 2, 
                            ifelse(control_d$final_Batch == "GWAS4", 1, NA))
control_d$GWAS1_5 <- ifelse(control_d$final_Batch == "GWAS1", 2, 
                            ifelse(control_d$final_Batch == "GWAS5", 1, NA))
control_d$GWAS1_GEM <- ifelse(control_d$final_Batch == "GWAS1", 2, 
                              ifelse(control_d$final_Batch == "GEM", 1, NA))

control_d$GWAS2_3 <- ifelse(control_d$final_Batch == "GWAS2", 2, 
                            ifelse(control_d$final_Batch == "GWAS3", 1, NA))
control_d$GWAS2_4 <- ifelse(control_d$final_Batch == "GWAS2", 2, 
                            ifelse(control_d$final_Batch == "GWAS4", 1, NA))
control_d$GWAS2_5 <- ifelse(control_d$final_Batch == "GWAS2", 2, 
                            ifelse(control_d$final_Batch == "GWAS5", 1, NA))
control_d$GWAS2_GEM <- ifelse(control_d$final_Batch == "GWAS2", 2, 
                              ifelse(control_d$final_Batch == "GEM", 1, NA))

control_d$GWAS3_4 <- ifelse(control_d$final_Batch == "GWAS3", 2, 
                            ifelse(control_d$final_Batch == "GWAS4", 1, NA))
control_d$GWAS3_5 <- ifelse(control_d$final_Batch == "GWAS3", 2, 
                            ifelse(control_d$final_Batch == "GWAS5", 1, NA))
control_d$GWAS3_GEM <- ifelse(control_d$final_Batch == "GWAS3", 2, 
                              ifelse(control_d$final_Batch == "GEM", 1, NA))

control_d$GWAS4_5 <- ifelse(control_d$final_Batch == "GWAS4", 2, 
                            ifelse(control_d$final_Batch == "GWAS5", 1, NA))
control_d$GWAS4_GEM <- ifelse(control_d$final_Batch == "GWAS4", 2, 
                              ifelse(control_d$final_Batch == "GEM", 1, NA))

control_d$GWAS5_GEM <- ifelse(control_d$final_Batch == "GWAS5", 2, 
                              ifelse(control_d$final_Batch == "GEM", 1, NA))

pseudoGWAS_pheno <- control_d[,c("V1","V1","GWAS1_2","GWAS1_3","GWAS1_4","GWAS1_5","GWAS1_GEM",
                                 "GWAS2_3","GWAS2_4","GWAS2_5","GWAS2_GEM","GWAS3_4","GWAS3_5","GWAS3_GEM",
                                 "GWAS4_5","GWAS4_GEM","GWAS5_GEM")]
pseudoGWAS_covar <- control_d[,c("V1","V1","Sex","Race")]

colnames(pseudoGWAS_pheno) <- c("FID","IID", 
                                "GWAS1_2","GWAS1_3","GWAS1_4","GWAS1_5","GWAS1_GEM",
                                "GWAS2_3","GWAS2_4","GWAS2_5","GWAS2_GEM","GWAS3_4","GWAS3_5","GWAS3_GEM",
                                "GWAS4_5","GWAS4_GEM","GWAS5_GEM")
colnames(pseudoGWAS_covar) <- c("FID","IID","Sex","Race")

table(pseudoGWAS_covar$Sex)
sum(is.na(pseudoGWAS_covar$Sex))

write.table(pseudoGWAS_pheno, "./pseudoGWAS_input/AllCohorts_Controls_pseudoGWAS_pheno.txt", 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(pseudoGWAS_covar, "./pseudoGWAS_input/AllCohorts_Controls_pseudoGWAS_covar.txt", 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

###
pseudoGWAS_pheno <- read.table("./pseudoGWAS_input/AllCohorts_Controls_pseudoGWAS_pheno.txt", header = TRUE)
pseudoGWAS_covar <- read.table("./pseudoGWAS_input/AllCohorts_Controls_pseudoGWAS_covar.txt", header = TRUE)

# recode sex and race 
pseudoGWAS_covar$sex_num_recode <- ifelse(pseudoGWAS_covar$Sex == "F", 2, 1)
# table(pseudoGWAS_covar$Race)
pseudoGWAS_covar$race_num_recode <- ifelse(pseudoGWAS_covar$Race %in% c("W", "White"), 1, 
                                           ifelse(pseudoGWAS_covar$Race %in% c("B", "Black"), 2, 
                                                  ifelse(pseudoGWAS_covar$Race == "A", 3, 
                                                         ifelse(pseudoGWAS_covar$Race < 0, NA, 4))))
table(pseudoGWAS_covar$race_num_recode)
sum(is.na(pseudoGWAS_covar$race_num_recode)) 
# 241

control_PCA <- fread("./pseudoGWAS_input/All_Cohorts_allChr_noDups_Controls_PCA.eigenvec")

# update pheno and covar - remove duplicates in plink before
pseudoGWAS_pheno <- pseudoGWAS_pheno[pseudoGWAS_pheno$IID %in% control_PCA$IID,]
pseudoGWAS_covar <- pseudoGWAS_covar[pseudoGWAS_covar$IID %in% control_PCA$IID,]

control_phenoPC <- merge(pseudoGWAS_covar, control_PCA, by = intersect(colnames(pseudoGWAS_covar), colnames(control_PCA)))

# table(control_phenoPC$race_num_recode)
control_racePC <- control_phenoPC[,c("FID", "IID", "race_num_recode", paste0("PC",1:20))]
cols <- rep(NA, nrow(control_phenoPC))

cols[control_racePC$race_num_recode == 1] <- "blue" # white
cols[control_racePC$race_num_recode == 2] <- "red" # black
cols[control_racePC$race_num_recode == 3] <- "green" # asian
cols[control_racePC$race_num_recode == 4] <- "orange" # others
cols[is.na(control_racePC$race_num_recode)] <- "grey" # NA

control_par.coord <- control_racePC[,c(paste0("PC",1:20))]
control_rangel <- apply(control_par.coord, 2, function(x) range(x)[1])
control_rangeh <- apply(control_par.coord, 2, function(x) range(x)[2])
control_std.coord <- control_par.coord
for (i in 1:14)
  control_std.coord[,i] <- (control_par.coord[,i] - control_rangel[i])/(control_rangeh[i]-control_rangel[i])

pdf("./pseudoGWAS_input/All_Cohorts_Control_N7555_PCA Parallel Coorinates.pdf")
plot(c(0,15), c(0,1), type = 'n', axes = FALSE, ylab = "", xlab = "",
     main = "Parallel Coordinates Plot for GWAS1,2,3,4,5,GEM cohorts controls:
     White: blue, Black: red, Asian: green, Others: orange")
for (j in 1:13)
  for (i in sample(1:nrow(control_std.coord)) )
    lines(c(j,j+1), control_std.coord[i,c(j,j+1)], col=cols[i], lwd=0.25)
axis(1, at = 1:14, labels = paste("PC",1:14, sep = "."))
dev.off()

# from the plot, use PC8

control_covar_d <- control_phenoPC[,c("FID", "IID", "sex_num_recode", paste0("PC",1:8))]
str(control_covar_d)
write.table(control_covar_d, "./pseudoGWAS_input/All_Cohorts_Control_Covariates_for_PLINK.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(pseudoGWAS_pheno, "./pseudoGWAS_input/All_Cohorts_Control_Phenotypes_for_PLINK.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# save.image("All_Cohort_Case_Control_Prep.RData")
load("All_Cohort_Case_Control_Prep.RData")
