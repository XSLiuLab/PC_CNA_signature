# Load packages -----------------------------------------------------------

library(sigminer)
library(tidyverse)
library(NMF)


# Reading data ------------------------------------------------------------
Maf = data.table::fread("/public/data/maf/all.maf")
# Remove all NA columns
Maf = Maf[,which(unlist(lapply(Maf, function(x)!all(is.na(x))))),with=F]
Maf = read_maf(Maf)
save(Maf, file="output/PRAD_TCGA_plus_dbGap_Maf.RData")

# Prepare data and estimate number -------------------------------------------
ncores = 12

system.time(
  Maf.derive <- sig_derive(Maf, cores = ncores, ref_genome = "BSgenome.Hsapiens.UCSC.hg38",
                        useSyn = TRUE)
)
save(Maf.derive, file = "output/PRAD_TCGA_plus_dbGap_Maf.derive.RData")

load(file = "output/PRAD_TCGA_plus_dbGap_Maf.derive.RData")
# Remove the effect of hyper mutated samples (not removing hyper-mutated samples)
nmf_matrix = handle_hyper_mutation(Maf.derive$nmf_matrix)

# EST.Maf = sig_estimate(Maf.derive$nmf_matrix,
#                        range = 2:10, nrun = 50, cores = ncores, use_random = TRUE,
#                        save_plots = FALSE,
#                        verbose = TRUE)
# save(EST.Maf, file = "output/EST.PRAD_TCGA_plus_dbGap_Maf.RData")

EST.Maf.rm_hyper = sig_estimate(nmf_matrix,
                                range = 2:10, nrun = 50, cores = ncores, use_random = TRUE,
                                save_plots = FALSE,
                                verbose = TRUE)

save(EST.Maf.rm_hyper, file = "output/EST.PRAD_TCGA_plus_dbGap_Maf_rm_hyper.RData")
load(file = "output/EST.PRAD_TCGA_plus_dbGap_Maf_rm_hyper.RData")

show_sig_number_survey(EST.Maf)
show_sig_number_survey2(EST.Maf$survey, EST.Maf$survey.random)


show_sig_number_survey(EST.Maf.rm_hyper)
show_sig_number_survey2(EST.Maf.rm_hyper$survey, EST.Maf.rm_hyper$survey.random, what = "sparseness")
EST.Maf.rm_hyper$survey_plot


# Extract signatures ------------------------------------------------------

Sig.SNV = sig_extract(nmf_matrix, n_sig = 6, nrun = 50, cores = ncores)
save(Sig.SNV, file = "output/Sig.PRAD_TCGA_plus_dbGap_rm_hyper.RData")

show_sig_profile(Sig.SNV, mode = "mutation")
get_sig_similarity(Sig.SNV)
show_cosmic_sig_profile(c(3, 26, 7, 30, 15, 1))

dd = get_groups(Sig.SNV)
table(dd$enrich_sig)

show_sig_exposure(Sig.SNV, rm_space = T, cutoff = 2000)
