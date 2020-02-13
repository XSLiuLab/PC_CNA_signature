# NOTE: the raw data are big, thus not maintained in repository
library(dplyr)
library(stringr)
source("code/99-functions.R")

# Obtaining absolute copy number from Sequenza ------------------------------

extract_seqz_cnv("/home/wsx/projects/prad_signature/cnv_calling/sequenza/seqz_wes_result",
                 target_path = "data/CNV_from_sequenza.tsv"
)

# Test
seqz_cnv <- readr::read_tsv("data/CNV_from_sequenza.tsv")
test <- sigminer::read_copynumber(seqz_cnv,
                                  genome_build = "hg38",
                                  complement = FALSE, verbose = TRUE
)

# Obtaining purity and ploidy ---------------------------------------------

extract_seqz_purity_and_ploidy("/home/wsx/projects/prad_signature/cnv_calling/sequenza/seqz_wes_result",
                               target_path = "data/PRAD_Purity_and_Ploidy_Sequenza.tsv"
)

# Transform data to GISTIC input
# and store as .seg file
seqz_to_GISTIC2("/home/wsx/projects/prad_signature/cnv_calling/sequenza/seqz_wes_result",
                  target_path = 'data/PRAD_CNA_hg38_seqz.seg',
                  rm_samps = "WCMC160-SRR3146971")

