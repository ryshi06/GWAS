
library(topr)
library(data.table)
library(qqman)
library(dplyr)
library(Haplin)

setwd("/ix/kfan/Ruyu/multi-trait_GWAS/")

###--------------------------------------------------------------------------###
### Set 1 Abeta
###--------------------------------------------------------------------------###

# on GRCh38
### read in multi-trait gwas file 
# gwas1 <- fread("./output/ABCDS_multitrait_Set1_AB_residuals.assoc.txt")
gwas1 <- fread("./output/ABCDS_multitrait_Set1_AB_residuals_MAF01.assoc.txt")

head(gwas1)
gwas1$uniqID <- paste0(gwas1$chr, "_", gwas1$ps)
gwas1$direction_1 <- ifelse(gwas1$beta_1 > 0, "+", "-")
gwas1$direction_2 <- ifelse(gwas1$beta_2 > 0, "+", "-")
gwas1$direction_3 <- ifelse(gwas1$beta_3 > 0, "+", "-")
gwas1$direction_4 <- ifelse(gwas1$beta_4 > 0, "+", "-")
gwas1$direction <- paste(gwas1$direction_1, gwas1$direction_2, gwas1$direction_3, gwas1$direction_4, sep="")

# filter by MAF 5%
gwas1_MAF05 <- gwas1[gwas1$af >= 0.05,] %>% 
  select(chr, rs, ps, p_lrt, direction) %>%
  rename(CHROM = chr, rsID = rs, POS = ps, P = p_lrt)

# calculate lambda
gwas1_MAF05$CHISQ <- qchisq(gwas1_MAF05$P,1,lower.tail=FALSE)

# annotation
gwas1_MAF05_sig <- gwas1_MAF05[gwas1_MAF05$P < 1e-6, ] %>% annotate_with_nearest_gene()
gwas1_MAF05_sig <- merge(gwas1_MAF05_sig, gwas1[,c("rs", "allele1", "allele0", "af")], by.x = "rsID", by.y = "rs", all.x = TRUE)

lambda <- median(qchisq(1-gwas1_MAF05$P, 1))/qchisq(0.5,1)
m <- length(gwas1_MAF05$P)
x <- rev(qchisq(1-(1:m)/(m), 1))
y <- gwas1_MAF05$CHISQ

png("./Figures/set1_AB_MAF05_GEMMA_QQ_N199.png",width=480,height=480)
pQQ(gwas1_MAF05$P, mark=F, cex.lab=1.25, main="Set1-AB: \nAmyloid-PET Centiloid, Plasma Ab40, Plasma Ab42, Ab42/40 Ratio")
mtext(paste0("lamda = ",round(lambda,3)), side=3, line=0)
dev.off()

# rename for Manhattan plot
gwas1_MAF05 <- gwas1_MAF05 %>% rename(CHR = CHROM, BP = POS, SNP = rsID)
png("./Figures/set1_AB_GEMMA_MAF05_Manhattan.png", width = 1080, height = 480)
par(mar=c(4,5,3,1))
manhattan(gwas1_MAF05, col=c("black","green2","dodgerblue4","darkslategray1","darkviolet","red","yellow","gray"), 
          cex=0.7, cex.axis=1.5, cex.lab=1.5, genomewideline=-log10(5e-08), ylim=c(0,10), 
          suggestiveline = -log10(1e-6), main="Set1-AB: \nAmyloid-PET Centiloid, Plasma Ab40, Plasma Ab42, Ab42/40 Ratio")
dev.off()

### get statistics of shared significant SNPs

# read in individual gwas file
#Ab40
x <- fread("./Results/ABC_DS_GWAS/ABCDS_Amy_Ab40.assoc.linear.frq.gene.nonemissing", header=TRUE)
x1 <- subset(x, select=c("CHR", "BP", "SNP", "P", "MAF", "BETA", "GENE"))

#Ab42
x <- fread("./Results/ABC_DS_GWAS/ABCDS_Amy_Ab42.assoc.linear.frq.gene.nonemissing", header=TRUE)
x2 <- subset(x, select=c("CHR", "BP", "SNP", "P", "MAF", "BETA", "GENE"))

#Ab4042 Ratio
x <- fread("./Results/ABC_DS_GWAS/ABCDS_Norm_Ab4240_ratio.assoc.linear.frq.gene.nonemissing", header=TRUE)
x3 <- subset(x, select=c("CHR", "BP", "SNP", "P", "MAF", "BETA", "GENE"))

#Centiloid
x <- fread("./Results/ABC_DS_GWAS/ABCDS_Centiloid.assoc.linear.frq.gene.nonemissing", header=TRUE)
x4 <- subset(x, select=c("CHR", "BP", "SNP", "P", "MAF", "BETA", "GENE"))

AB_shared_snps <- c(unique(x1[x1$P < 5e-08, ]$SNP),
                    unique(x2[x2$P < 5e-08, ]$SNP),
                    unique(x3[x3$P < 5e-08, ]$SNP),
                    unique(x4[x4$P < 5e-08, ]$SNP),
                    unique(gwas1_MAF05[gwas1_MAF05$P < 5e-08,]$SNP))
length(AB_shared_snps)

ab_sig <- list()
for (file in c(paste0("x",1:4))){
  x <- get(file)
  ab_sig[[file]] <- x[x$SNP %in% AB_shared_snps, c("SNP", "P", "MAF","BETA")]
}

ab_sig1 <- do.call(cbind, ab_sig)

colnames(ab_sig1) <- gsub("x1.", "Ab40.", colnames(ab_sig1))
colnames(ab_sig1) <- gsub("x2.", "Ab42.", colnames(ab_sig1))
colnames(ab_sig1) <- gsub("x3.", "Ab4240.", colnames(ab_sig1))
colnames(ab_sig1) <- gsub("x4.", "Centiloid.", colnames(ab_sig1))

ab_sig1 <- ab_sig1[,c("Ab40.SNP", "Ab40.MAF", "Ab40.BETA", "Ab40.P", "Ab42.MAF", "Ab42.BETA", "Ab42.P",
                      "Ab4240.MAF", "Ab4240.BETA", "Ab4240.P","Centiloid.MAF", "Centiloid.BETA", "Centiloid.P")]
colnames(ab_sig1) <- c("rs", "Ab40.MAF", "Ab40.BETA", "Ab40.P", "Ab42.MAF", "Ab42.BETA", "Ab42.P",
                       "Ab4240.MAF", "Ab4240.BETA", "Ab4240.P","Centiloid.MAF", "Centiloid.BETA", "Centiloid.P")

ab_sig2 <- merge(ab_sig1, 
                 gwas1[gwas1$rs %in% AB_shared_snps, c("chr","ps","rs","allele1","allele0","af","p_lrt","direction")],
                 by = "rs")
dim(ab_sig2)
# [1]  8 20

ab_sig3 <- ab_sig2[order(ab_sig2$chr, ab_sig2$ps), 
                   c("rs", "chr", "ps", "allele1", "allele0", "af", "p_lrt", "direction", 
                     "Ab40.MAF", "Ab40.BETA", "Ab40.P","Ab42.MAF", "Ab42.BETA", "Ab42.P",
                     "Ab4240.MAF", "Ab4240.BETA", "Ab4240.P", "Centiloid.MAF", "Centiloid.BETA", "Centiloid.P")] %>%
  dplyr::rename(CHROM = chr, SNP = rs, POS = ps, P = p_lrt)

ab_sig_annotated <- ab_sig3 %>% annotate_with_nearest_gene()

write.table(ab_sig_annotated, "./Results/AB_5E-08_SNPs_annotated.txt", sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)

##ab_sig_annotated %>% 
##  gt::gt() %>%
##  gt::tab_header(
##    title = "Abeta Genome-wide Significant SNPs",
##    subtitle = "Multi-trait GWAS + single-trait GWAS"
##  ) %>%
##  gt::fmt_number(
##    columns = -c(SNP, CHROM, POS, allele1, allele0, direction, P, Ab42.P, Ab40.P, Ab4240.P, Centiloid.P),
##    decimals = 2,
##    use_seps = FALSE
##  ) %>%
##  gt::fmt_scientific(
##    columns = c(P, Ab42.P, Ab40.P, Ab4240.P, Centiloid.P),
##    decimals = 2
##  ) %>%
##  gt::tab_spanner(label = "Ab40", columns = starts_with("Ab40.")) %>%
##  gt::tab_spanner(label = "Ab42", columns = starts_with("Ab42.")) %>%
##  gt::tab_spanner(label = "Ab4240", columns = starts_with("Ab4240.")) %>%
##  gt::tab_spanner(label = "Centiloid", columns = starts_with("Centiloid.")) %>%
##  gt::cols_label(
##    Ab40.MAF = "MAF", Ab40.BETA = "BETA", Ab40.P = "P",
##    Ab42.MAF = "MAF", Ab42.BETA = "BETA", Ab42.P = "P",
##    Ab4240.MAF = "MAF", Ab4240.BETA = "BETA", Ab4240.P = "P",
##    Centiloid.MAF = "MAF", Centiloid.BETA = "BETA", Centiloid.P = "P"
##  ) %>%
##  gt::cols_width(c(SNP, biotype) ~ px(160),
##                 c(ends_with(".P"), P, POS, Gene_Symbol) ~ px(120),
##                 everything() ~ px(70)) %>%
##  gt::gtsave(filename = "./Figures/AB_5E-08_SNPs_annotated.png", 
##             vwidth = 3000, vheight = 1800, expand = 4)

###--------------------------------------------------------------------------###
### Set 2 Tau
###--------------------------------------------------------------------------###

gwas2 <- fread("./output/ABCDS_multitrait_Set2_Tau_residuals_MAF01.assoc.txt")

head(gwas2)
gwas2$uniqID <- paste0(gwas2$chr, "_", gwas2$ps)
gwas2$direction_1 <- ifelse(gwas2$beta_1 > 0, "+", "-")
gwas2$direction_2 <- ifelse(gwas2$beta_2 > 0, "+", "-")
gwas2$direction_3 <- ifelse(gwas2$beta_3 > 0, "+", "-")
gwas2$direction_4 <- ifelse(gwas2$beta_4 > 0, "+", "-")
gwas2$direction <- paste(gwas2$direction_1, gwas2$direction_2, gwas2$direction_3, gwas2$direction_4, sep="")

# filter by MAF 5%
gwas2_MAF05 <- gwas2[gwas2$af >= 0.05,] %>% 
  select(chr, rs, ps, p_lrt, direction) %>%
  dplyr::rename(CHROM = chr, rsID = rs, POS = ps, P = p_lrt)

# calculate lambda
gwas2_MAF05$CHISQ <- qchisq(gwas2_MAF05$P,1,lower.tail=FALSE)

# annotation
gwas2_MAF05_sig <- gwas2_MAF05[gwas2_MAF05$P < 1e-06, ] %>% annotate_with_nearest_gene()
gwas2_MAF05_sig <- merge(gwas2_MAF05_sig, gwas2[,c("rs", "allele1", "allele0", "af")], by.x = "rsID", by.y = "rs", all.x = TRUE)

lambda <- median(qchisq(1-gwas2_MAF05$P, 1))/qchisq(0.5,1)
m <- length(gwas2_MAF05$P)
x <- rev(qchisq(1-(1:m)/(m), 1))
y <- gwas2_MAF05$CHISQ

png("./Figures/set2_Tau_MAF05_GEMMA_QQ_N76.png",width=480,height=480)
pQQ(gwas1_MAF05$P, mark=F, cex.lab=1.25, main="Set2-Tau: \nTau-PET, Plasma Tau181, Plasma Tau217, Plasma Total Tau")
mtext(paste0("lamda = ",round(lambda,3)), side=3, line=0)
dev.off()

png("./Figures/set2_Tau_GEMMA_MAF05_Manhattan.png", width = 1080, height = 600)
manhattan(gwas2_MAF05, col=c("black","green2","dodgerblue4","darkslategray1","darkviolet","red","yellow","gray"), 
          cex=0.7, cex.axis=1.5, cex.lab=1.5, genomewideline=-log10(5e-08), ylim=c(0,10), 
          suggestiveline = -log10(1e-6), main="Set2-Tau (N = 76): \nTau-PET, Plasma Tau181, Plasma Tau217, Plasma Total Tau")
dev.off()

###--------------------------------------------------------------------------###
### Set 3 Tau (with imputed value)
###--------------------------------------------------------------------------###

##gwas3 <- fread("./output/ABCDS_multitrait_Set3_Tau_residuals_impu.assoc.txt")

##head(gwas3)
##gwas3$uniqID <- paste0(gwas3$chr, "_", gwas3$ps)
##gwas3$direction_1 <- ifelse(gwas3$beta_1 > 0, "+", "-")
##gwas3$direction_2 <- ifelse(gwas3$beta_2 > 0, "+", "-")
##gwas3$direction_3 <- ifelse(gwas3$beta_3 > 0, "+", "-")
##gwas3$direction_4 <- ifelse(gwas3$beta_4 > 0, "+", "-")
##gwas3$direction <- paste(gwas3$direction_1, gwas3$direction_2, gwas3$direction_3, gwas3$direction_4, sep="")

# filter by MAF 5%
##gwas3_MAF05 <- gwas3[gwas3$af >= 0.05,] %>% 
##  select(chr, rs, ps, p_lrt, direction) %>%
##  rename(CHROM = chr, rsID = rs, POS = ps, P = p_lrt)

# calculate lambda
##gwas3_MAF05$CHISQ <- qchisq(gwas3_MAF05$P,1,lower.tail=FALSE)

# annotation
##gwas3_MAF05_sig <- gwas3_MAF05[gwas3_MAF05$P < 1e-06, ] %>% annotate_with_nearest_gene()
##gwas3_MAF05_sig <- merge(gwas3_MAF05_sig, gwas3[,c("rs", "allele1", "allele0", "af")], by.x = "rsID", by.y = "rs", all.x = TRUE)

##lambda <- median(qchisq(1-gwas3_MAF05$P, 1))/qchisq(0.5,1)
##m <- length(gwas3_MAF05$P)
##x <- rev(qchisq(1-(1:m)/(m), 1))
##y <- gwas3_MAF05$CHISQ

##png("./Figures/set3_Tau_imputed_MAF05_GEMMA_QQ_N320.png")
##pQQ(gwas3_MAF05$P, mark=F, cex.lab=1.25, main="Set3-Tau (imputed): \nTau-PET, Plasma Tau181, Plasma Tau217, Plasma Total Tau")
##mtext(paste0("lamda = ",round(lambda,3)), side=3, line=0)
##dev.off()

# rename for Manhattan plot
##gwas3_MAF05 <- gwas3_MAF05 %>% rename(CHR = CHROM, BP = POS, SNP = rsID)

##png("./Figures/set2_set3_Tau_GEMMA_MAF05_Manhattan.png", width = 1080, height = 600)
##par(mfrow=c(2,1))
##par(mar=c(0,5,3,3))
##manhattan(gwas2_MAF05, col=c("black","green2","dodgerblue4","darkslategray1","darkviolet","red","yellow","gray"), 
##          cex=0.7, cex.axis=1.5, cex.lab=1.5, genomewideline=-log10(5e-08), ylim=c(0,10), 
##          suggestiveline = -log10(1e-6), main="Set2-Tau: \nTau-PET, Plasma Tau181, Plasma Tau217, Plasma Total Tau")
##par(mar=c(5,5,3,3))
##manhattan(gwas3_MAF05, col=c("black","green2","dodgerblue4","darkslategray1","darkviolet","red","yellow","gray"), 
##          cex=0.7, cex.axis=1.5, cex.lab=1.5, genomewideline=-log10(5e-08), ylim=c(10,0), xlab="", xaxt="n", 
##          suggestiveline = -log10(1e-6))
##          # main="Set3-Tau (imputed): \nTau-PET, Plasma Tau181, Plasma Tau217, Plasma Total Tau")
##dev.off()

### pick the SNPs passed the suggestive significant line in Set2 and Set3
##shared_snps <- c(unique(gwas2_MAF05_sig$rsID), unique(gwas3_MAF05_sig$rsID))
##shared_sig <- merge(gwas2[gwas2$rs %in% shared_snps, c("chr","ps","rs","allele1","allele0","af","p_lrt","direction")], 
##                    gwas3[gwas3$rs %in% shared_snps, c("rs","af","p_lrt","direction")], by = "rs")
##colnames(shared_sig) <- gsub(".x",".tau", colnames(shared_sig))
##colnames(shared_sig) <- gsub(".y",".tau.imputed", colnames(shared_sig))
##write.table(shared_sig[order(shared_sig$chr, shared_sig$ps),], "./Results/Set2_Set3_Tau_1E06_overlapped_SNPs.txt", sep = "\t", 
##            col.names = TRUE, row.names = FALSE, quote = FALSE)

# read in individual gwas file
#inv_norm_Tau_PET
x <- fread("./Results/ABC_DS_GWAS/ABCDS_Norm_Tau.assoc.linear.frq.gene.nonemissing", header=TRUE)
x1 <- subset(x, select=c("CHR", "BP", "SNP", "P", "MAF", "BETA", "GENE")) %>%
  filter(MAF >= 0.01)

#inv_norm_Tau_181
x <- fread("./Results/ABC_DS_GWAS/ABCDS_Nor_pTau181.assoc.linear.frq.gene.nonemissing", header=TRUE)
x2 <- subset(x, select=c("CHR", "BP", "SNP", "P", "MAF", "BETA", "GENE")) %>%
  filter(MAF >= 0.01)

#inv_norm_Tau_217
x <- fread("./Results/ABC_DS_GWAS/ABCDS_Nor_pTau217.assoc.linear.frq.gene.nonemissing", header=TRUE)
x3 <- subset(x, select=c("CHR", "BP", "SNP", "P", "MAF", "BETA", "GENE")) %>%
  filter(MAF >= 0.01)

#inv_norm_Totl_Tau
x <- fread("./Results/ABC_DS_GWAS/ABCDS_Nor_TotalpTau.assoc.linear.frq.gene.nonemissing", header=TRUE)
x4 <- subset(x, select=c("CHR", "BP", "SNP", "P", "MAF", "BETA", "GENE")) %>%
  filter(MAF >= 0.01)

Tau_shared_snps <- c(unique(x1[x1$P < 5e-08, ]$SNP),
                     unique(x2[x2$P < 5e-08, ]$SNP),
                     unique(x3[x3$P < 5e-08, ]$SNP),
                     unique(x4[x4$P < 5e-08, ]$SNP),
                     unique(gwas2_MAF05[gwas2_MAF05$P < 1e-06,]$rsID))
length(Tau_shared_snps)

# found no overlaps of Tau_shared_snps in individual GWAS and multi-trait GWAS
tau_sig <- list()
for (file in c(paste0("x",1:4))){
  x <- get(file)
  tau_sig[[file]] <- x[x$SNP %in% unique(Tau_shared_snps), c("SNP", "P", "MAF","BETA")]
}

tau_sig1 <- do.call(cbind, tau_sig)

colnames(tau_sig1) <- gsub("x1.", "Tau_PET.", colnames(tau_sig1))
colnames(tau_sig1) <- gsub("x2.", "Tau_181.", colnames(tau_sig1))
colnames(tau_sig1) <- gsub("x3.", "Tau_217.", colnames(tau_sig1))
colnames(tau_sig1) <- gsub("x4.", "Totl_Tau.", colnames(tau_sig1))

tau_sig1 <- tau_sig1[,c("Tau_PET.SNP", "Tau_PET.MAF", "Tau_PET.BETA", "Tau_PET.P", "Tau_181.MAF", "Tau_181.BETA", "Tau_181.P",
                        "Tau_217.MAF", "Tau_217.BETA", "Tau_217.P","Totl_Tau.MAF", "Totl_Tau.BETA", "Totl_Tau.P")]
colnames(tau_sig1) <- c("rs", "Tau_PET.MAF", "Tau_PET.BETA", "Tau_PET.P", "Tau_181.MAF", "Tau_181.BETA", "Tau_181.P",
                        "Tau_217.MAF", "Tau_217.BETA", "Tau_217.P","Totl_Tau.MAF", "Totl_Tau.BETA", "Totl_Tau.P")
# head(shared_sig)
tau_sig2 <- merge(tau_sig1, gwas2, by = "rs", all.x = TRUE)
# tau_sig2_sub <- tau_sig2[tau_sig2$P < 1e-06,]
dim(tau_sig2)
# [1]  9 42

tau_sig3 <- tau_sig2[order(tau_sig2$chr, tau_sig2$ps), 
                   c("rs", "chr", "ps", "allele1", "allele0", "af", "p_lrt", "direction", 
                     "Tau_PET.MAF", "Tau_PET.BETA", "Tau_PET.P", "Tau_181.MAF", "Tau_181.BETA", "Tau_181.P",
                     "Tau_217.MAF", "Tau_217.BETA", "Tau_217.P","Totl_Tau.MAF", "Totl_Tau.BETA", "Totl_Tau.P")] %>%
  dplyr::rename(CHROM = chr, SNP = rs, POS = ps, P = p_lrt) 

tau_sig_annotated <- tau_sig3 %>% annotate_with_nearest_gene()

write.table(tau_sig_annotated, "./Results/Tau_1E-06_SNPs_annotated.txt", sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)

##tau_sig_annotated %>% 
##  filter(!is.na(P)) %>% 
##  gt::gt() %>%
##  gt::tab_header(
##    title = "Tau Genome-wide Suggestive Significant SNPs",
##    subtitle = "Multi-trait GWAS + single-trait GWAS"
##  ) %>%
##  gt::fmt_number(
##    columns = -c(SNP, CHROM, POS, allele1, allele0, direction, P),
##    decimals = 2,
##    use_seps = FALSE
##  ) %>%
##  gt::fmt_scientific(
##    columns = P,
##    decimals = 2
##  ) %>%
##  gt::tab_spanner(label = "Tau_PET", columns = starts_with("Tau_PET.")) %>%
##  gt::tab_spanner(label = "Tau_181", columns = starts_with("Tau_181.")) %>%
##  gt::tab_spanner(label = "Tau_217", columns = starts_with("Tau_217.")) %>%
##  gt::tab_spanner(label = "Totl_Tau", columns = starts_with("Totl_Tau.")) %>%
##  gt::cols_label(
##    Tau_PET.MAF = "MAF", Tau_PET.BETA = "BETA", Tau_PET.P = "P",
##    Tau_181.MAF = "MAF", Tau_181.BETA = "BETA", Tau_181.P = "P",
##    Tau_217.MAF = "MAF", Tau_217.BETA = "BETA", Tau_217.P = "P",
##    Totl_Tau.MAF = "MAF", Totl_Tau.BETA = "BETA", Totl_Tau.P = "P"
##  ) %>%
##  gt::cols_width(c(SNP, biotype) ~ px(170),
##                 c(ends_with(".P"), P, POS, Gene_Symbol) ~ px(120),
##                 everything() ~ px(70)) %>%
##  gt::gtsave(filename = "./Figures/Tau_1E-06_SNPs_annotated.png", 
##             vwidth = 3000, vheight = 1800, expand = 4)

