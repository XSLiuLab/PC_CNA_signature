# Combine and clean IDs deeply
# using IDs from:
# 1) FACETS result
# 2) 02-tidy_clinical_data.R
# 3) clinical data from Nat.Gen 2018 (1013 samples)
library(tidyverse)


CNV = read_copynumber("data/CNV_from_dbGAP_PLUS_TCGA_WES_CVAL150.tsv", genome_build = "hg38",
                      complement = FALSE, verbose = TRUE)

data_clincal <- readxl::read_excel("../data/NIHMS938028-supplement-5.xlsx", skip = 2)
names(data_clincal)[2] <- "Tumor_Sample_Barcode"
