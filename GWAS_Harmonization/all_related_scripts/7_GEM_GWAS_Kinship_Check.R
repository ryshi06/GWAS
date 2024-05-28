
###--------------------------------------------------------------------------###
###
### GEM GWAS QC check
###
###--------------------------------------------------------------------------###

rm(list=ls())
setwd("/home/Data/GEM_GWAS/Frank")
options(stringsAsFactors = FALSE)

###--------------------------------------------------------------------------###

ibd <- read.csv("GEM_GWAS_Kinship_Estimate.txt", sep="\t")
load("GEM_GWAS_Annotation.RData")
d <- pData(scanAnnot)

# str(ibd)
# str(d)

ibd$rel.dup <- ifelse(ibd$exp.rel=="Dup" & ibd$obs.rel=="Dup", 1, 0)
ibd2 <- subset(ibd, ibd$rel.dup==0)

# str(ibd2)

###--------------------------------------------------------------------------###

ids <- unique(sort(c(ibd2$Individ1, ibd2$Individ2)))

d2 <- subset(d, corrected_scanName%in%ids, 
             select=c("corrected_scanName", "plate", "well_position", 
                      "sex", "race", "age", "center"))

x <- do.call("rbind",lapply(ibd2$Individ1, function(x){
  d2[which(d2$corrected_scanName==x),]}))
  
# str(x)
# str(ibd2)
# identical(x$corrected_scanName,as.character(ibd2$Individ1))

names(x) <- c("Ind1","plate_Ind1","well_pos_Ind1","sex_Ind1",
              "race_Ind1", "age_Ind1", "center_Ind1")

ibd3 <- cbind(ibd2,x)
# str(ibd3)

y <- do.call("rbind",lapply(ibd2$Individ2, function(x){
  d2[which(d2$corrected_scanName==x),]}))

# str(y)
# str(ibd3)
# identical(y$corrected_scanName,as.character(ibd3$Individ2))

names(y) <- c("Ind2","plate_Ind2","well_pos_Ind2","sex_Ind2",
              "race_Ind2", "age_Ind2", "center_Ind2")

ibd4 <- cbind(ibd3,y)
# str(ibd4)

# identical(as.character(ibd4$Individ1), ibd4$Ind1)
# identical(as.character(ibd4$Individ2), ibd4$Ind2)

ibdf <- subset(ibd4, select=c("ii", "exp.rel", "obs.rel", "Ind1", "age_Ind1", "sex_Ind1",
                              "race_Ind1", "center_Ind1", "plate_Ind1", "well_pos_Ind1",
                              "Ind2", "age_Ind2", "sex_Ind2", "race_Ind2", "center_Ind2",
                              "plate_Ind2", "well_pos_Ind2"))

write.table(ibdf, "GEM_GWAS_Kinship_Check.txt", sep="\t", quote = FALSE, row.names = FALSE)

###--------------------------------------------------------------------------###














