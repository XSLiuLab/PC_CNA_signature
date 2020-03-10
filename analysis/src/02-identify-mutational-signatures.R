# Load packages -----------------------------------------------------------

library(sigminer)
library(tidyverse)
library(NMF)


# Reading data ------------------------------------------------------------
Maf <- data.table::fread("/public/data/maf/all.maf")
# Remove all NA columns
Maf <- Maf[, which(unlist(lapply(Maf, function(x) !all(is.na(x))))), with = F]
# openxlsx::write.xlsx(Maf, file = "output/PRAD_MAF.xlsx")

Maf <- read_maf(Maf)
save(Maf, file = "output/PRAD_TCGA_plus_dbGap_Maf.RData")

load(file = "output/PRAD_TCGA_plus_dbGap_Maf.RData")
# Prepare data and estimate number -------------------------------------------
ncores <- 12

system.time(
  Maf.tally <- sig_tally(Maf,
    cores = ncores, ref_genome = "BSgenome.Hsapiens.UCSC.hg38",
    useSyn = TRUE
  )
)
save(Maf.tally, file = "output/PRAD_TCGA_plus_dbGap_Maf.tally.RData")

load(file = "output/PRAD_TCGA_plus_dbGap_Maf.tally.RData")
# Remove the effect of hyper mutated samples (not removing hyper-mutated samples)
nmf_matrix <- handle_hyper_mutation(Maf.tally$nmf_matrix)

EST.Maf <- sig_estimate(nmf_matrix,
  range = 2:10, nrun = 50, cores = ncores, use_random = TRUE,
  save_plots = FALSE,
  verbose = TRUE
)

save(EST.Maf, file = "output/EST.PRAD_TCGA_plus_dbGap_Maf.RData")
save(nmf_matrix, file = "output/Maf_matrix.RData")

load(file = "output/EST.PRAD_TCGA_plus_dbGap_Maf.RData")

show_sig_number_survey(EST.Maf)
show_sig_number_survey2(EST.Maf$survey, EST.Maf$survey.random)

# Extract signatures ------------------------------------------------------
Sig.SNV <- sig_extract(nmf_matrix, n_sig = 3, nrun = 50, cores = ncores)
save(Sig.SNV, file = "output/Sig.PRAD_TCGA_plus_dbGap_Maf.RData")

get_sig_similarity(Sig.SNV)
get_sig_similarity(Sig.SNV, sig_db = "SBS")

# Extract 5 signatures to compare with copy number signatures
Sig.SNV5 <- sig_extract(nmf_matrix, n_sig = 5, nrun = 50, cores = 20)
save(Sig.SNV5, file = "output/Sig5.Maf.RData")

# # Auto converge to 3
# Sig.SNV.auto = sig_auto_extract(nmf_matrix, nrun = 100, destdir = "output/BayesianNMF_MutSig",
#                                 cores = 10, recover = TRUE)
