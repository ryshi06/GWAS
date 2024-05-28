# cohort get overlap SNPs 

rm(list=ls())
options(stringsAsFactors = FALSE)

# set the working directory to ***your own*** directory 
setwd("/zfs1/User/cohort/")

# laod the packages that we are going to use later
# This is the main package we used in GWAS QC
library(GWASTools)
# This package can load Excel (.xlsx) files into R 
library(gdata)

# a small function to select which "x" is not in "y"
# which will be using lot of times later
"%nin%" <- function(x,y){
    !x%in%y
}

###############################################################
### reference file obtained from https://www.strand.org.uk/ ###
###############################################################

### replace the actual .strand file in the following command line
ref <- read.csv("/zfs1/User/GWAS_Chip_Info/HumanOmni2-5-8-v1-1-C-b37-strand/HumanOmni2-5-8-v1-1-C-b37.strand", 
                sep="\t", header = FALSE)

colnames(ref) <- c("Name", "Chr", "Position", "%match_to_genome", "strand", "TOP")

# select the columns we need, discard others
d1 <- subset(ref, select=c("Name", "Chr", "Position"))

# rename the columns 
names(d1) <- c("snpName", "chromosome", "position")

#############################################
### reference file obtained from Illumina ###
#############################################

### replace the actual _LocusReport.txt file in the following command line
# ref <- read.csv("/zfs1/User/GWAS_Chip_Info/GDA-8V1-0/infinium-global-diversity-array-8-v1-0_D1_LocusReport.txt", 
#                sep="\t")

# select the columns we need, discard others
# d1 <- subset(ref,select=c("Name", "Chr", "Position"))

# rename the columns 
# names(d1) <- c("snpName", "chromosome", "position")

### make the chromosome index as an integer: (following statement is from "DataCleaning.pdf") 
#   chromosome, an integer mapping for each chromosome, with values 1-27, mapped in order
#   from 1-22, 23=X, 24=XY (the pseudoautosomal region), 25=Y, 26=M (the mitochondrial probes), 
#   and 27=U (probes with unknown positions).                                                                      caution. See the manual pages for more details.)

d1$chromosome[d1$chromosome=="X"] <- 23
d1$chromosome[d1$chromosome=="Y"] <- 25
d1$chromosome[d1$chromosome=="MT"] <- 26
d1$chromosome[d1$chromosome=="0"] <- 27

# If the gene annotation file ("Multi-EthnicGlobal_D1.annotated.txt") is going to be used, 
# then use 23 for XY instead of 24
# d1$chromosome[d1$chromosome=="XY"] <- 23
d1$chromosome[d1$chromosome=="XY"] <- 24

# make the chromosome as a interger
d1$chromosome <- as.integer(d1$chromosome)

############################################################
### get overlapped variants of reference and sample file ###
############################################################

### replace the actual .txt sample file in the following command line

# read one sample file
sample <- read.table("./Raw_Data/cohort_Sample_sampleName.txt", skip = 9, header = TRUE)

# get overlapped SNPs of two files
d1_overlap <- d1[d1$snpName %in% intersect(sample$snp, d1$snpName), ]
d1_overlap$Index <- 1:nrow(d1_overlap)

write.table(d1_overlap, "cohort_Ref_Overlap.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# loop over cohort Raw_Data/ to generate sample files with overlap SNPs
cohort_file_list <- list.files("./Raw_Data", pattern = "^cohort_Sample_")
dir.create("./Processed_Overlap", showWarnings = FALSE)

start_time <- Sys.time()
count = 1
for (file in cohort_file_list) {
    
    # print the file name
    print(paste(count, file))
    # read the file
    data <- read.table(file.path("./Raw_Data/", file), skip = 9, header = TRUE)
    data_overlap <- data[data$snp %in% d1_overlap$snpName,]
    
    # define the output file path
    output_file <- file.path("./Processed_Overlap", file)
    
    # save the filtered data to the subfolder
    write.table(data_overlap, output_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    
    count = count + 1
    
    # once loop 100 times, run garbage collection to release memory
    if (count %% 100 == 0) {
        gc()
    }
    
}
end_time <- Sys.time()
print(paste("Elapsed time:", end_time-start_time))
