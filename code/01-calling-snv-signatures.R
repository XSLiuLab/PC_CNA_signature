library(sigminer)
library(NMF)
library(data.table)


# Reading data ------------------------------------------------------------
Maf = data.table::fread("/public/data/maf/all.maf")
# Remove all NA columns
Maf = Maf[,which(unlist(lapply(Maf, function(x)!all(is.na(x))))),with=F]
Maf = read_maf(Maf)
save(Maf, file="data/PRAD_Maf.RData")

# Prepare data and estimate rank -------------------------------------------
ncores = 12

system.time(
  Maf.Pre <- sig_derive(Maf, cores = ncores, ref_genome = "BSgenome.Hsapiens.UCSC.hg38",
                    useSyn = TRUE)
)

save(Maf.Pre, file = "output/Maf.Pre.RData")

# Auto-extract
Sig.Bayesian.Maf = sig_auto_extract(Maf.Pre$nmf_matrix, result_prefix = "BayesNMF_Maf", nrun = 100,
                                     destdir = "output/signature", cores = 16)
# Sig.Bayesian.Maf = sig_auto_extract(Maf.Pre$nmf_matrix, result_prefix = "BayesNMF_Maf", nrun = 100,
#                                     destdir = "output/signature", cores = 16, recover = T)

# Extract by hand
Maf_Est = sig_estimate(Maf.Pre$nmf_matrix,
                            range = 2:10, nrun = 50, cores = ncores, use_random = TRUE,
                            save_plots = FALSE,
                            verbose = TRUE)

save(Maf_Est, file = "output/Maf_Est.RData")

show_rank_survey(Maf_Est)
# Select best signature number
Sig.SNV= sig_extract(Maf.Pre$nmf_matrix, n_sig = 3, nrun = 100, cores = ncores)
saveRDS(Sig.SNV, file = "output/NMF_snv_signature.rds")

# Normalise by Row
show_sig_profile(Sig.SNV, mode = "mutation")

Sig.SNV = readRDS("output/NMF_snv_signature.rds")
get_sig_similarity(Sig.SNV)
get_sig_similarity(Sig.SNV, sig_db = "SBS")

