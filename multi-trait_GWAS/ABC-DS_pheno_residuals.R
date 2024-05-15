#################################
### ABC-DS MultiPhen Analysis ###
#################################

###-----------------------------------------------###
### set up working directory and load R libraries ###
###-----------------------------------------------###

rm(list=ls())
setwd("/ix/kfan/Ruyu/multi-trait_GWAS")
options(stringsAsFactors = FALSE)

library(data.table)
# gc()

### ------------------------------------------------------------------------ ###
### ABC-DS phenotype preparation
###
### Dec 20th, 2023
### ------------------------------------------------------------------------ ###

pheno_dir <- "./Data/Phenotype/"
geno_dir <- "./Data/Genotype/"

output_dir <- "./Results/"

### ------------------------------------------------------------------------ ###
### Functions
### ------------------------------------------------------------------------ ###

inverse.normalize <- function(i) {
  stopifnot(is.vector(i, mode = "numeric"))
  qnorm((rank(i,
              na.last = "keep",
              ties.method = "random"
  ) - 0.5) / sum(!is.na(i)))
}

count_missing_values <- function(data) {
  result_df <- data.frame(Column = names(data), Missing_Count = numeric(length(data)))
  
  for (i in 1:ncol(data)) {
    result_df$Missing_Count[i] <- sum(is.na(data[, i]))
  }
  return(result_df)
}

### ------------------------------------------------------------------------ ###
# read pheno file

pheno <- read.table(paste0(pheno_dir, "ABC-DS_Geno_Pheno.txt"), sep = "\t", header = TRUE)

# str(pheno)
colnames(pheno)

### ------------------------------------------------------------------------ ###
### get residuals of required phenotype sets
### ------------------------------------------------------------------------ ###

phenos <- c("centiloid_value","Ab40", "Ab42", "inv_norm_AB4240",
           "inv_norm_Tau_PET","inv_norm_Tau_181","inv_norm_Tau_217","inv_norm_Totl_Tau")

covar <- c("age", "Sex", "base_dementia","PC1","PC2","PC3","PC4")

### ------------------------------------------------------------------------ ###
### check missing patterns
### ------------------------------------------------------------------------ ###

# check missing patterns
missing <- list()

pheno_dd <- pheno[,phenos]
covar_dd <- pheno[,covar]

for (phe in phenos){
  missing[[phe]] <- ifelse(is.na(pheno_dd[, phe]),1,0)
}

missing_dd <- do.call(cbind, missing)

# run logistic regression to check MAR assumption
covar_others <- covar_dd[,3:5]
p_vals <- list()

for (phe in phenos){
  for (col in 1:ncol(covar_others)){
    tmp  <- cbind(missing_dd[,phe], covar_others[,col])
    logistic_reg <- summary(glm(tmp[,1] ~ tmp[,2]), 
                            family = binomial)
    p_vals[[phe]][col] <- logistic_reg$coefficients[2,4]
  }
}

### ------------------------------------------------------------------------ ###
## set 1: AB
### ------------------------------------------------------------------------ ###

set1_pheno <- c("centiloid_value","Ab40", "Ab42", "inv_norm_AB4240")
set1 <- subset(pheno, !is.na(IID) & !is.na(centiloid_value) & !is.na(Ab40) & !is.na(Ab42) & !is.na(inv_norm_AB4240) & !is.na(base_dementia) & !is.na(age) & !is.na(Sex),
               select = c("FID", "IID", covar, set1_pheno))
# dim(set1)
# [1] 199  13

# table(set1$base_dementia)
#     0     1 
#   163   36

# count_missing_values(set1)

# set1_residuals 
set1_residuals <- cbind(set1[,c("FID", "IID")], as.data.frame(matrix(nrow=nrow(set1), ncol=length(set1_pheno))))

for (i in 1:length(set1_pheno)){
  model <- lm(set1[, which(colnames(set1) == set1_pheno[i])] ~ 
                age + Sex + base_dementia + PC1 + PC2 + PC3 + PC4, data = set1)
  print(set1_pheno[i])
  set1_residuals[,(i+2)] <- residuals(model)
}
colnames(set1_residuals) <- c("FID", "IID", set1_pheno)
# head(set1_residuals)

set1_residuals <- subset(set1_residuals, !is.na(IID) & (!is.na(centiloid_value) | !is.na(Ab40) | !is.na(Ab42) | !is.na(inv_norm_AB4240)), )

## for (residual in 1:4){
##   print(shapiro.test(set1_residuals[,2+residual]))
## }

set1_phenos <- set1[set1$IID %in% set1_residuals$IID, c("FID", "IID", set1_pheno)]

## for (pheno in 1:4){
##   print(shapiro.test(set1_phenos[,2+pheno]))
## }

set1_covars <- set1[set1$IID %in% set1_residuals$IID, c("age", "Sex", "base_dementia","PC1","PC2","PC3","PC4")]

write.table(set1_residuals[,c("FID", "IID")], paste0(pheno_dir, "ABCDS_multitrait_Set1_AB_FID_IID.txt"), 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(set1_residuals, paste0(pheno_dir, "ABCDS_multitrait_Set1_AB_residuals.txt"), 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(set1_phenos, paste0(pheno_dir, "ABCDS_multitrait_Set1_AB_phenos.txt"), 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(set1_covars, paste0(pheno_dir, "ABCDS_multitrait_Set1_AB_covars.txt"), 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

### ------------------------------------------------------------------------ ###
## set 2: Tau
### ------------------------------------------------------------------------ ###

# "Tau_PET", "IN_Tau_PET", "pTau181", "IN_Tau_181", "pTau217", "IN_Tau_217", "TotalpTau", "IN_Total_Tau"
set2_pheno <- c("inv_norm_Tau_PET","inv_norm_Tau_181","inv_norm_Tau_217","inv_norm_Totl_Tau")
set2 <- subset(pheno, !is.na(IID) & !is.na(inv_norm_Tau_PET) & !is.na(inv_norm_Tau_181) & !is.na(inv_norm_Tau_217) & !is.na(inv_norm_Totl_Tau) & !is.na(base_dementia) & !is.na(age) & !is.na(Sex),
               select = c("FID", "IID", covar, set2_pheno))

# dim(set2)
# [1] 76 14

# table(set2$base_dementia)
#  0   1 
# 69   7

# count_missing_values(set2)

###--- set2_residuals 
set2_residuals <- cbind(set2[,c("FID", "IID")], as.data.frame(matrix(nrow=length(set2$IID), ncol=length(set2_pheno))))

for (i in 1:length(set2_pheno)){
  model <- lm(set2[, which(colnames(set2) == set2_pheno[i])] ~ 
                age + Sex + base_dementia + PC1 + PC2 + PC3 + PC4, data = set2)
  print(set2_pheno[i])
  set2_residuals[,(i+2)] <- residuals(model)
}
colnames(set2_residuals) <- c("FID", "IID", set2_pheno)
# head(set2_residuals)

set2_residuals <- subset(set2_residuals, !is.na(IID) & (!is.na(inv_norm_Tau_PET) | !is.na(inv_norm_Tau_181) | !is.na(inv_norm_Tau_217) | !is.na(inv_norm_Totl_Tau)), )

set2_phenos <- set2[set2$IID %in% set2_residuals$IID, c("FID", "IID", set2_pheno)]
set2_covars <- set2[set2$IID %in% set2_residuals$IID, c("age", "Sex", "base_dementia","PC1","PC2","PC3","PC4")]

write.table(set2_residuals[,c("FID", "IID")], paste0(pheno_dir, "ABCDS_multitrait_Set2_Tau_FID_IID.txt"), 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(set2_residuals, paste0(pheno_dir, "ABCDS_multitrait_Set2_Tau_residuals.txt"), 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(set2_phenos, paste0(pheno_dir, "ABCDS_multitrait_Set2_Tau_phenos.txt"), 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(set2_covars, paste0(pheno_dir, "ABCDS_multitrait_Set2_Tau_covars.txt"), 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

### ------------------------------------------------------------------------ ###
### set 3: Tau with missing
### ------------------------------------------------------------------------ ###

# "Tau_PET", "IN_Tau_PET", "pTau181", "IN_Tau_181", "pTau217", "IN_Tau_217", "TotalpTau", "IN_Total_Tau"
set3_pheno <- c("inv_norm_Tau_PET","inv_norm_Tau_181","inv_norm_Tau_217","inv_norm_Totl_Tau")
set3 <- subset(pheno, !is.na(IID) & (!is.na(inv_norm_Tau_PET) | !is.na(inv_norm_Tau_181) | !is.na(inv_norm_Tau_217) | !is.na(inv_norm_Totl_Tau)), select = c("FID", "IID", covar, set3_pheno))

rownames(set3) <- 1:nrow(set3)

# dim(set3)
# [1] 333  13

# table(set3$base_dementia)
#   0   1 
# 231  89 
# sum(is.na(set3$base_dementia))
# 13

# count_missing_values(set3)

###--- set3_residuals 
set3_residuals <- cbind(set3[,c("FID", "IID")], as.data.frame(matrix(nrow=length(set3$IID), ncol=length(set3_pheno))))

for (i in 1:length(set3_pheno)){
  model <- lm(set3[, which(colnames(set3) == set3_pheno[i])] ~ 
                age + Sex + base_dementia + PC1 + PC2 + PC3 + PC4, data = set3)
  print(set3_pheno[i])
  
  res <- as.data.frame(residuals(model))
  res$row_number <- rownames(res)
  colnames(res) <- c("res", "row_number")
  missing_rows <- setdiff(1:nrow(set3), res$row_number)
  
  if(length(missing_rows) == 0){
    set3_residuals[,(i+2)] <- residuals(model)
  } else {
    res <- as.data.frame(residuals(model))
    res$row_number <- rownames(res)
    colnames(res) <- c("res", "row_number")
    missing_rows <- setdiff(1:nrow(set3), res$row_number)
    # assign NA to those missing rows
    res2 <- rbind(res, data.frame(res = NA, row_number = missing_rows))
    # reorder the rows
    res2 <- res2[order(as.numeric(res2$row_number)), ]
    set3_residuals[,(i+2)] <- res2$res
  }
}

colnames(set3_residuals) <- c("FID", "IID", set3_pheno)
# head(set3_residuals)
# dim(set3_residuals)
set3_residuals <- subset(set3_residuals, !is.na(IID) & (!is.na(inv_norm_Tau_PET) | !is.na(inv_norm_Tau_181) | !is.na(inv_norm_Tau_217) | !is.na(inv_norm_Totl_Tau)), )

set3_phenos <- set3[set3$IID %in% set3_residuals$IID, c("FID", "IID", set3_pheno)]
set3_covars <- set3[set3$IID %in% set3_residuals$IID, c("age", "Sex", "base_dementia","PC1","PC2","PC3","PC4")]

write.table(set3_residuals[,c("FID", "IID")], paste0(pheno_dir, "ABCDS_multitrait_Set3_Tau_FID_IID.txt"), 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(set3_residuals, paste0(pheno_dir, "ABCDS_multitrait_Set3_Tau_residuals.txt"), 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(set3_phenos, paste0(pheno_dir, "ABCDS_multitrait_Set3_Tau_phenos.txt"), 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(set3_covars, paste0(pheno_dir, "ABCDS_multitrait_Set3_Tau_covars.txt"), 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

### ------------------------------------------------------------------------ ###
### set 4: Abeta (no PET)
### ------------------------------------------------------------------------ ###

# "Tau_PET", "IN_Tau_PET", "pTau181", "IN_Tau_181", "pTau217", "IN_Tau_217", "TotalpTau", "IN_Total_Tau"
set4_pheno <- c("Ab40", "Ab42", "inv_norm_AB4240")
set4 <- subset(pheno, !is.na(IID) & !is.na(Ab40) | !is.na(Ab42) | !is.na(inv_norm_AB4240), select = c("FID", "IID", covar, set4_pheno))
set4 <- set4[!is.na(set4$base_dementia),]

rownames(set4) <- 1:nrow(set4)

# dim(set4)
# [1] 275 12

# table(set4$base_dementia)
#   0   1 
# 204  71 
# sum(is.na(set4$base_dementia))

# count_missing_values(set4)

###--- set4_residuals 
set4_residuals <- cbind(set4[,c("FID", "IID")], as.data.frame(matrix(nrow=length(set4$IID), ncol=length(set4_pheno))))

for (i in 1:length(set4_pheno)){
  model <- lm(set4[, which(colnames(set4) == set4_pheno[i])] ~ 
                age + Sex + base_dementia + PC1 + PC2 + PC3 + PC4, data = set4)
  print(set4_pheno[i])
  
  res <- as.data.frame(residuals(model))
  res$row_number <- rownames(res)
  colnames(res) <- c("res", "row_number")
  missing_rows <- setdiff(1:nrow(set4), res$row_number)
  
  if(length(missing_rows) == 0){
    set4_residuals[,(i+2)] <- residuals(model)
  } else {
    res <- as.data.frame(residuals(model))
    res$row_number <- rownames(res)
    colnames(res) <- c("res", "row_number")
    missing_rows <- setdiff(1:nrow(set4), res$row_number)
    # assign NA to those missing rows
    res2 <- rbind(res, data.frame(res = NA, row_number = missing_rows))
    # reorder the rows
    res2 <- res2[order(as.numeric(res2$row_number)), ]
    set4_residuals[,(i+2)] <- res2$res
  }
}

colnames(set4_residuals) <- c("FID", "IID", set4_pheno)
# head(set4_residuals)

set4_residuals <- subset(set4_residuals, !is.na(IID) & !is.na(Ab40) | !is.na(Ab42) | !is.na(inv_norm_AB4240), )

set4_phenos <- set4[set4$IID %in% set4_residuals$IID, c("FID", "IID", set4_pheno)]
set4_covars <- set4[set4$IID %in% set4_residuals$IID, c("age", "Sex", "base_dementia","PC1","PC2","PC3","PC4")]

write.table(set4_residuals[,c("FID", "IID")], paste0(pheno_dir, "ABCDS_multitrait_Set4_AB_noPET_FID_IID.txt"), 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(set4_residuals, paste0(pheno_dir, "ABCDS_multitrait_Set4_AB_noPET_residuals.txt"), 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(set4_phenos, paste0(pheno_dir, "ABCDS_multitrait_Set4_AB_noPET_phenos.txt"), 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(set4_covars, paste0(pheno_dir, "ABCDS_multitrait_Set4_AB_noPET_covars.txt"), 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

### ------------------------------------------------------------------------ ###
## set 5: Tau + PET
### ------------------------------------------------------------------------ ###

# "centiloid_value", "Tau_PET", "IN_Tau_PET", "pTau181", "IN_Tau_181", "pTau217", "IN_Tau_217", "TotalpTau", "IN_Total_Tau"
set5_pheno <- c("centiloid_value", "inv_norm_Tau_PET","inv_norm_Tau_181","inv_norm_Tau_217","inv_norm_Totl_Tau")
set5 <- subset(pheno, !is.na(IID) & !is.na(centiloid_value) | !is.na(inv_norm_Tau_PET) | !is.na(inv_norm_Tau_181) | !is.na(inv_norm_Tau_217) | !is.na(inv_norm_Totl_Tau), select = c("FID", "IID", covar, set5_pheno))

rownames(set5) <- 1:nrow(set5)

# dim(set5)
# [1] 337  14

# table(set5$base_dementia)
#   0   1 
# 234  90 
# sum(is.na(set5$base_dementia))
# 13

# count_missing_values(set5)

###--- set5_residuals 
set5_residuals <- cbind(set5[,c("FID", "IID")], as.data.frame(matrix(nrow=length(set5$IID), ncol=length(set5_pheno))))

for (i in 1:length(set5_pheno)){
  model <- lm(set5[, which(colnames(set5) == set5_pheno[i])] ~ 
                age + Sex + base_dementia + PC1 + PC2 + PC3 + PC4, data = set5)
  print(set5_pheno[i])
  
  res <- as.data.frame(residuals(model))
  res$row_number <- rownames(res)
  colnames(res) <- c("res", "row_number")
  missing_rows <- setdiff(1:nrow(set5), res$row_number)
  
  if(length(missing_rows) == 0){
    set5_residuals[,(i+2)] <- residuals(model)
  } else {
    res <- as.data.frame(residuals(model))
    res$row_number <- rownames(res)
    colnames(res) <- c("res", "row_number")
    missing_rows <- setdiff(1:nrow(set5), res$row_number)
    # assign NA to those missing rows
    res2 <- rbind(res, data.frame(res = NA, row_number = missing_rows))
    # reorder the rows
    res2 <- res2[order(as.numeric(res2$row_number)), ]
    set5_residuals[,(i+2)] <- res2$res
  }
}

colnames(set5_residuals) <- c("FID", "IID", set5_pheno)
# head(set5_residuals)

# set5_residuals <- subset(set5_residuals, !(is.na(centiloid_value) & is.na(inv_norm_Tau_PET) & is.na(inv_norm_Tau_181) & is.na(inv_norm_Tau_217) & is.na(inv_norm_Totl_Tau)))
dim(set5_residuals)
# [1] 337  7

set5_phenos <- set5[set5$IID %in% set5_residuals$IID, c("FID", "IID", set5_pheno)]
set5_covars <- set5[set5$IID %in% set5_residuals$IID, c("age", "Sex", "base_dementia","PC1","PC2","PC3","PC4")]

write.table(set5_residuals[,c("FID", "IID")], paste0(pheno_dir, "ABCDS_multitrait_Set5_Tau_plus_PET_FID_IID.txt"), 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(set5_residuals, paste0(pheno_dir, "ABCDS_multitrait_Set5_Tau_plus_PET_residuals.txt"), 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(set5_phenos, paste0(pheno_dir, "ABCDS_multitrait_Set5_Tau_plus_PET_phenos.txt"), 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(set5_covars, paste0(pheno_dir, "ABCDS_multitrait_Set5_Tau_plus_PET_covars.txt"), 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

### ------------------------------------------------------------------------ ###
## set 6: Tau + PET (no missing)
### ------------------------------------------------------------------------ ###

# "centiloid_value", "Tau_PET", "IN_Tau_PET", "pTau181", "IN_Tau_181", "pTau217", "IN_Tau_217", "TotalpTau", "IN_Total_Tau"
set6_pheno <- c("centiloid_value", "inv_norm_Tau_PET","inv_norm_Tau_181","inv_norm_Tau_217","inv_norm_Totl_Tau")
set6 <- subset(pheno, !is.na(IID) & !is.na(centiloid_value) | !is.na(inv_norm_Tau_PET) | !is.na(inv_norm_Tau_181) | !is.na(inv_norm_Tau_217) | !is.na(inv_norm_Totl_Tau), select = c("FID", "IID", covar, set6_pheno))

rownames(set6) <- 1:nrow(set6)

# dim(set6)
# [1] 337  14

# table(set6$base_dementia)
#   0   1 
# 234  90 
# sum(is.na(set6$base_dementia))
# 13

# count_missing_values(set6)

###--- set6_residuals 
set6_residuals <- cbind(set6[,c("FID", "IID")], as.data.frame(matrix(nrow=length(set6$IID), ncol=length(set6_pheno))))

for (i in 1:length(set6_pheno)){
  model <- lm(set6[, which(colnames(set6) == set6_pheno[i])] ~ 
                age + Sex + base_dementia + PC1 + PC2 + PC3 + PC4, data = set6)
  print(set6_pheno[i])
  
  res <- as.data.frame(residuals(model))
  res$row_number <- rownames(res)
  colnames(res) <- c("res", "row_number")
  missing_rows <- setdiff(1:nrow(set6), res$row_number)
  
  if(length(missing_rows) == 0){
    set6_residuals[,(i+2)] <- residuals(model)
  } else {
    res <- as.data.frame(residuals(model))
    res$row_number <- rownames(res)
    colnames(res) <- c("res", "row_number")
    missing_rows <- setdiff(1:nrow(set6), res$row_number)
    # assign NA to those missing rows
    res2 <- rbind(res, data.frame(res = NA, row_number = missing_rows))
    # reorder the rows
    res2 <- res2[order(as.numeric(res2$row_number)), ]
    set6_residuals[,(i+2)] <- res2$res
  }
}

colnames(set6_residuals) <- c("FID", "IID", set6_pheno)
# head(set6_residuals)

set6_residuals <- subset(set6_residuals, 
                         !is.na(IID) & (!is.na(centiloid_value) & !is.na(inv_norm_Tau_PET) & !is.na(inv_norm_Tau_181) & !is.na(inv_norm_Tau_217) & !is.na(inv_norm_Totl_Tau)), )
dim(set6_residuals)
# [1] 71  7

set6_phenos <- set6[set6$IID %in% set6_residuals$IID, c("FID", "IID", set6_pheno)]
set6_covars <- set6[set6$IID %in% set6_residuals$IID, c("age", "Sex", "base_dementia","PC1","PC2","PC3","PC4")]

write.table(set6_residuals[,c("FID", "IID")], paste0(pheno_dir, "ABCDS_multitrait_set6_Tau_plus_PET_FID_IID.txt"), 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(set6_residuals, paste0(pheno_dir, "ABCDS_multitrait_set6_Tau_plus_PET_residuals.txt"), 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(set6_phenos, paste0(pheno_dir, "ABCDS_multitrait_set6_Tau_plus_PET_phenos.txt"), 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(set6_covars, paste0(pheno_dir, "ABCDS_multitrait_set6_Tau_plus_PET_covars.txt"), 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# save.image(paste0(pheno_dir, "ABCDS_multitrait_Pheno.RData"))
load(paste0(pheno_dir, "ABCDS_multitrait_Pheno.RData"))

################################################################################
################################################################################
### update .fam file by adding phenotypes/residuals
################################################################################
################################################################################

setwd("/ix/kfan/Ruyu/multi-trait_GWAS")

library(data.table)

### ------------------------------------------------------------------------ ###
### Set 1: Abeta group (complete) (N = 199)
### ------------------------------------------------------------------------ ###

fam <- read.table("./Data/Genotype/ABCDS_multitrait_Set1_AB.fam")
table(fam$V5)
head(fam)
residuals <- fread("./Data/Phenotype/ABCDS_multitrait_Set1_AB_residuals.txt")
residuals.fam <- merge(fam, residuals, by = c("V1","V2"))
write.table(residuals.fam[,which(colnames(residuals.fam) != "V6.x")], "./Data/Genotype/ABCDS_multitrait_Set1_AB_residuals.fam", 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

phenos <- fread("./Data/Phenotype/ABCDS_multitrait_Set1_AB_phenos.txt")
phenos.fam <- merge(fam, phenos, by = c("V1","V2"))
write.table(phenos.fam[,which(colnames(phenos.fam) != "V6.x")], "./Data/Genotype/ABCDS_multitrait_Set1_AB_phenos.fam", 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

### ------------------------------------------------------------------------ ###
### Set 2: Tau group (complete) (N = 76)
### ------------------------------------------------------------------------ ###

fam <- read.table("./Data/Genotype/ABCDS_multitrait_Set2_Tau.fam")
table(fam$V5)
head(fam)
residuals <- fread("./Data/Phenotype/ABCDS_multitrait_Set2_Tau_residuals.txt")
residuals.fam <- merge(fam, residuals, by = c("V1","V2"))
write.table(residuals.fam[,which(colnames(residuals.fam) != "V6.x")], "./Data/Genotype/ABCDS_multitrait_Set2_Tau_residuals.fam", 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

phenos <- fread("./Data/Phenotype/ABCDS_multitrait_Set2_Tau_phenos.txt")
phenos.fam <- merge(fam, phenos, by = c("V1","V2"))
write.table(phenos.fam[,which(colnames(phenos.fam) != "V6.x")], "./Data/Genotype/ABCDS_multitrait_Set2_Tau_phenos.fam", 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

### ------------------------------------------------------------------------ ###
### Set 3: Tau group (N = 320)
### ------------------------------------------------------------------------ ###

fam <- read.table("./Data/Genotype/ABCDS_multitrait_Set3_Tau.fam")
table(fam$V5)
head(fam)
residuals <- fread("./Data/Phenotype/ABCDS_multitrait_Set3_Tau_residuals.txt")
residuals.fam <- merge(fam, residuals, by = c("V1","V2"))
write.table(residuals.fam[,which(colnames(residuals.fam) != "V6.x")], "./Data/Genotype/ABCDS_multitrait_Set3_Tau_residuals.fam", 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

phenos <- fread("./Data/Phenotype/ABCDS_multitrait_Set3_Tau_phenos.txt")
phenos.fam <- merge(fam, phenos, by = c("V1","V2"))
write.table(phenos.fam[,which(colnames(phenos.fam) != "V6.x")], "./Data/Genotype/ABCDS_multitrait_Set3_Tau_phenos.fam", 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

### update .fam again with predicted residuals 
## fam <- read.table("./Data/Genotype/ABCDS_multitrait_Set3_Tau.fam")
## table(fam$V5)
## head(fam)
## residuals <- fread("./output/ABCDS_multitrait_Set3_Tau_residuals_prdt.prdt.txt")
## colnames(residuals) <- paste0("R",1:5)
## residuals.fam <- cbind(fam, residuals)
## colnames(residuals.fam)
## write.table(residuals.fam[,which(!(colnames(residuals.fam) %in% c("V6", "R5")))], "./Data/Genotype/ABCDS_multitrait_Set3_Tau_residuals_prdt.fam", 
##             col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

### ------------------------------------------------------------------------ ###
### Set 4: Abeta (no PET) N = 275
### ------------------------------------------------------------------------ ###

fam <- read.table("./Data/Genotype/ABCDS_multitrait_Set4_AB_noPET.fam")
table(fam$V5)
head(fam)
residuals <- fread("./Data/Phenotype/ABCDS_multitrait_Set4_AB_noPET_residuals.txt")
residuals.fam <- merge(fam, residuals, by = c("V1","V2"))
write.table(residuals.fam[,which(colnames(residuals.fam) != "V6")], "./Data/Genotype/ABCDS_multitrait_Set4_AB_noPET_residuals.fam", 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

phenos <- fread("./Data/Phenotype/ABCDS_multitrait_Set4_AB_noPET_phenos.txt")
phenos.fam <- merge(fam, phenos, by = c("V1","V2"))
write.table(phenos.fam[,which(colnames(phenos.fam) != "V6")], "./Data/Genotype/ABCDS_multitrait_Set4_AB_noPET_phenos.fam", 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

### ------------------------------------------------------------------------ ###
### Set 5: Tau (plus PET) N = 326
### ------------------------------------------------------------------------ ###

fam <- read.table("./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET.fam")
table(fam$V5)
head(fam)
residuals <- fread("./Data/Phenotype/ABCDS_multitrait_Set5_Tau_plus_PET_residuals.txt")
residuals.fam <- merge(fam, residuals, by = c("V1","V2"))
write.table(residuals.fam[,which(colnames(residuals.fam) != "V6.x")], "./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET_residuals.fam", 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

phenos <- fread("./Data/Phenotype/ABCDS_multitrait_Set5_Tau_plus_PET_phenos.txt")
phenos.fam <- merge(fam, phenos, by = c("V1","V2"))
write.table(phenos.fam[,which(colnames(phenos.fam) != "V6.x")], "./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET_phenos.fam", 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

### update .fam again with predicted residuals 
## fam <- read.table("./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET.fam")
## table(fam$V5)
## head(fam)
## residuals <- fread("./output/ABCDS_multitrait_Set5_Tau_plus_PET_residuals_prdt.prdt.txt")
## colnames(residuals) <- paste0("R",1:6)
## residuals.fam <- cbind(fam, residuals)
## colnames(residuals.fam)
## write.table(residuals.fam[,which(!(colnames(residuals.fam) %in% c("V6", "R6")))], "./Data/Genotype/ABCDS_multitrait_Set5_Tau_plus_PET_residuals_prdt.fam", 
##             col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

### ------------------------------------------------------------------------ ###
### Set 6: Tau (plus PET) N = 71
### ------------------------------------------------------------------------ ###

fam <- read.table("./Data/Genotype/ABCDS_multitrait_set6_Tau_plus_PET.fam")
table(fam$V5)
head(fam)
residuals <- fread("./Data/Phenotype/ABCDS_multitrait_set6_Tau_plus_PET_residuals.txt")
residuals.fam <- merge(fam, residuals, by = c("V1","V2"))
write.table(residuals.fam[,which(colnames(residuals.fam) != "V6.x")], "./Data/Genotype/ABCDS_multitrait_set6_Tau_plus_PET_residuals.fam", 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

phenos <- fread("./Data/Phenotype/ABCDS_multitrait_set6_Tau_plus_PET_phenos.txt")
phenos.fam <- merge(fam, phenos, by = c("V1","V2"))
write.table(phenos.fam[,which(colnames(phenos.fam) != "V6.x")], "./Data/Genotype/ABCDS_multitrait_set6_Tau_plus_PET_phenos.fam", 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")





