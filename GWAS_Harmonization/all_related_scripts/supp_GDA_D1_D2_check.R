D1 <- read.csv("/ihome/kfan/rus39/GDA-8v1-0_D1.csv", skip = 7, header = TRUE, fill = TRUE)
D2 <- read.csv("/ihome/kfan/rus39/GDA-8v1-0_D2.csv", skip = 7, header = TRUE, fill = TRUE)

all.equal(D1$IlmnID, D2$IlmnID)

D1[1:5,1:5]
# tail(D1[,1:5], n = 30)
D1_update <- D1[1:1904599,]
D2[1:5,1:5]
D2_update <- D2[1:1904599,]

comparison_result <- sapply(names(D1), function(col_name) {
  all.equal(D1[[col_name]], D2[[col_name]])
})

IGN235 <- read.table("/zfs1/kfan/GWAS4_Frank/GWAS4_Sample_IGN235.txt", skip = 9, header = TRUE, sep = "\t")
all.equal(D1_update$Name, IGN235$snp)



# check A ref and GWAS3 sample
HumanOmni_A <- read.table("/zfs1/kfan/Ruyu/harmonization_Sep5/GWAS_Chip_Info/HumanOmni2-5-8-v1-2-A-b37-strand/HumanOmni2.5-8-v1.2/HumanOmni2-5-8-v1-2-A-Loci-Report.txt",
                          sep = "\t", header = TRUE)
gwas3_s0359 <- read.table("/zfs1/kfan/Ruyu/harmonization_Sep5/GWAS3/Raw_Data/GWAS_3_Sample_s0359.txt",
                     sep = "\t", skip = 9, header = TRUE)

length(intersect(tolower(HumanOmni_A$Name), tolower(gwas3_s0359$snp)))
