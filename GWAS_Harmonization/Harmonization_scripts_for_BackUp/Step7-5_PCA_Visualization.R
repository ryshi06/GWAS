
library(ggplot2)
library(data.table)
library(dplyr)
library(openxlsx)

setwd("/ix/kfan/Ruyu/harmonization/All_Cohorts")

###--------------------------------------------------------------------------###
### load PCA output
###--------------------------------------------------------------------------###

eigenval <- fread("./QC_Check/AllCohorts_AllChr_PCA.eigenval")
eigenvec <- fread("./QC_Check/AllCohorts_AllChr_PCA.eigenvec")

###--------------------------------------------------------------------------###
### load phenotype
###--------------------------------------------------------------------------###

pheno <- read.xlsx("./Harmonization_Phenotype_2024_06_15.xlsx")
# length(intersect(pheno$scanName, eigenvec$IID))
table(pheno[is.na(pheno$Race),]$Batch)
# GWAS4 GWAS5 
#    15     1 

eigenvec_race <- merge(eigenvec, pheno[,c("scanName","Race")], by.x = "IID", by.y = "scanName") %>% distinct()
eigenvec_race$race_recode <- ifelse(is.na(eigenvec_race$Race), "Unknown", eigenvec_race$Race)
table(eigenvec_race$race_recode)
#  Asian   Black  Others Unknown   White 
#     25    1405      74      16    8911

##k <- 5  # number of clusters
##cluster_data <- eigenvec_race[, c("PC1", "PC2")]

##kmeans_result <- kmeans(cluster_data, centers = k)

# extract cluster assignments
##cluster_assignments <- kmeans_result$cluster

##eigenvec_race$race_recode <- factor(eigenvec_race$race_recode, levels = c("White", "Black", "Asian", "Others", "Unknown"))

##p <- ggplot(eigenvec_race) +
##        geom_point(aes(x = PC1, y = PC2, color = race_recode)) +
##        scale_color_manual(values = c("White"="#DA4453B2",
##                                      "Black"="#F6BA59B2",
##                                      "Asian"="#8BC163B2",
##                                      "Others"="#4B8AD6B2",
##                                      "Unknown"="#AAB2BCB2")) +
##        labs(title = "PCA plot of N = 10611") +
##        theme_bw()

### add a line to separate black and white
##library(e1071)
##eigenvec_bw <- subset(eigenvec_race, race_recode %in% c("White","Black"))
##rownames(eigenvec_bw) <- eigenvec_bw$IID
##eigenvec_bw <- eigenvec_bw[,c("PC1","PC2","race_recode")]

##eigenvec_bw$race_recode <- factor(eigenvec_bw$race_recode, levels = c("White","Black"))

# train SVM model
##svm_model <- svm(race_recode ~ ., data = eigenvec_bw, kernel = "linear", scale = FALSE, cost = 100000)

# create a grid of points to plot decision boundary
##x_grid <- seq(min(eigenvec_bw$PC1), max(eigenvec_bw$PC1), length.out = 100)
##y_grid <- seq(min(eigenvec_bw$PC2), max(eigenvec_bw$PC2), length.out = 100)
##grid <- expand.grid(PC1 = x_grid, PC2 = y_grid)

# predict class labels for grid points
##predicted_labels <- predict(svm_model, newdata = grid)

# get parameters of the hyperplane
##w <- t(svm_model$coefs) %*% svm_model$SV  # calculate the weight vector
##b <- -svm_model$rho

# get the slope/intercept of the linear line
### a <- -b/w[1,2]
### b <- -w[1,1]/w[1,2]

# plot the decision boundary
##plot(eigenvec_bw$PC1, eigenvec_bw$PC2, col = factor(eigenvec_bw$race_recode), pch = 19, main = "SVM Decision Boundary")
##points(grid$PC1, grid$PC2, col = factor(predicted_labels), pch = ".", cex = 1.5)
##legend("topright", legend = levels(eigenvec_bw$race_recode), col = c("White","Black"), pch = 19)

##plot(svm_model, data = eigenvec_bw)

##ggplot(eigenvec_bw) +
##geom_point(aes(x = PC1, y = PC2, color = race_recode)) +
##scale_color_manual(values = c("White"="#DA4453B2",
##"Black"="#F6BA59B2")) +
##theme_bw() +
##geom_abline(intercept = -b/w[1,2], slope = -w[1,1]/w[1,2], color = "black", linetype = "dashed")

# output the final figure
eigenvec_race$race_recode <- factor(eigenvec_race$race_recode, 
                                    levels = c("White", "Black", "Asian", "Others", "Unknown"))
pdf("./QC_Check/AllCohorts_AllChr_PCA.pdf", width = 10, height = 9)
ggplot(eigenvec_race) +
  geom_point(aes(x = PC1, y = PC2, color = race_recode, alpha = 0.8)) +
  scale_color_manual(values = c("White"="#DA4453B2",
                                "Black"="#F6BA59B2",
                                "Asian"="#8BC163B2",
                                "Others"="#4B8AD6B2",
                                "Unknown"="#AAB2BCB2")) +
  labs(title = "PCA plot of N = 10431") +
  theme_bw() +
  # geom_abline(intercept = -b/w[1,2], slope = -w[1,1]/w[1,2], color = "black", linetype = "dashed") +
  # geom_text_repel(data = extra_label[!(extra_label$IID %in% subset1$IID),],
  #                 aes(PC1, PC2, label = IID), color = "black", max.overlaps = 100) +
  # geom_text_repel(data = subset1[subset1$race_recode == "White",],
  #                 aes(PC1, PC2, label = IID), color = "darkred", max.overlaps = 100) +
  # geom_text_repel(data = subset1[subset1$race_recode == "Black",],
  #                 aes(PC1, PC2, label = IID), color = "darkorange", max.overlaps = 100) +
  geom_hline(yintercept = 0.05, color = "grey50", linetype = "dashed") +
  geom_vline(xintercept = 0.005, color = "grey50", linetype = "dashed")
dev.off()

###--------------------------------------------------------------------------###
### get the White/Black mis-classified sample lists
###--------------------------------------------------------------------------###

##eigenvec_race$Direction <- ifelse(eigenvec_race$IID %in% eigenvec_race[eigenvec_race$PC2 < -w[1,1]/w[1,2] * eigenvec_race$PC1 -b/w[1,2], ]$IID, "Down", "Up")
##table(eigenvec_race$Direction)
# Down   Up 
# 1433 9178 

##table(eigenvec_race$Direction, eigenvec_race$race_recode)

##nrow(eigenvec_race[eigenvec_race$Direction == "Down" & eigenvec_race$race_recode == "White",])
##nrow(eigenvec_race[eigenvec_race$Direction == "Up" & eigenvec_race$race_recode == "Black",])

white_inBlack <- eigenvec_race[eigenvec_race$PC1 > 0.005 & eigenvec_race$race_recode == "White",]
white_inBlack <- merge(white_inBlack[,c("IID", "PC1", "PC2")], pheno, by.x = "IID", by.y = "scanName")
dim(white_inBlack)
# [1] 25 16
white_inBlack$alert <- "self-reported White clustered in Black"
table(white_inBlack$Batch)
# GEM GWAS1 GWAS3 GWAS4 
#   4    12     4     5

black_inWhite <- eigenvec_race[eigenvec_race$PC1 < 0.005 & eigenvec_race$race_recode == "Black",]
black_inWhite <- merge(black_inWhite[,c("IID", "PC1", "PC2")], pheno, by.x = "IID", by.y = "scanName")
dim(black_inWhite)
# [1]  8 16
black_inWhite$alert <- "self-reported Black clustered in White"
table(black_inWhite$Batch)
# GWAS1 GWAS3 GWAS4 
#     1     1     6

asian_inWhite <- eigenvec_race[eigenvec_race$PC2 < 0.05 & eigenvec_race$race_recode == "Asian",]
asian_inWhite <- merge(asian_inWhite[,c("IID", "PC1", "PC2")], pheno, by.x = "IID", by.y = "scanName")
dim(asian_inWhite)
# [1]  6 16
asian_inWhite$alert <- "self-reported Asian clustered in White"
table(asian_inWhite$Batch)
# GWAS3 GWAS4 GWAS5 
#     1     4     1

problematic <- rbind(white_inBlack, black_inWhite, asian_inWhite) %>% distinct()
write.xlsx(problematic[order(problematic$Batch),], "./QC_Check/AllCohorts_AllChr_Harmonized_Race_Discrepancies.xlsx")

norace <- merge(eigenvec_race[,c("IID", "PC1", "PC2")], pheno[is.na(pheno$Race),], by.x = "IID", by.y = "scanName")
write.xlsx(norace[order(norace$Batch)], "./QC_Check/AllCohorts_AllChr_Harmonized_no_Race.xlsx")
