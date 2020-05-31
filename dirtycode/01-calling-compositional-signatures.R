library(tidyverse)
library(sigminer)
library(NMF)

load(file = "output/CNV.count.RData")
load(file = "output/Maf.Pre.RData")
Info <- readRDS("data/PRAD_CLINICAL.rds")

maps <- Info$CNV_ID
names(maps) <- Info$tumor_Run

CNV_mat <- CNV.count$nmf_matrix
SNV_mat <- Maf.Pre$nmf_matrix
rownames(SNV_mat)[startsWith(rownames(SNV_mat), "TCGA")] <- substr(rownames(SNV_mat)[startsWith(rownames(SNV_mat), "TCGA")], 1, 15)

rownames(SNV_mat) <- as.character(maps[rownames(SNV_mat)])

# Keep same IDs
same_ids <- intersect(rownames(CNV_mat), rownames(SNV_mat))
Mat <- cbind(CNV_mat[same_ids, ], SNV_mat[same_ids, ])


ncores <- 12
# Auto-extract
Sig.Bayesian.Comb <- sig_auto_extract(Mat,
  result_prefix = "BayesNMF_Comb", nrun = 100,
  destdir = "output/signature", cores = ncores
)


# Extract by hand
Comb_Est <- sig_estimate(Mat,
  range = 5:20, nrun = 100, cores = ncores, use_random = TRUE,
  save_plots = FALSE,
  verbose = TRUE, pConstant = 0.001
)

save(Comb_Est, file = "output/Comb_Est.RData")

show_rank_survey(Comb_Est)
# Select best signature number
Sig.Comb <- sig_extract(Mat, n_sig = 10, nrun = 100, cores = ncores, pConstant = 0.001)
saveRDS(Sig.Comb, file = "output/NMF_comb_signature.rds")

Sig.Comb <- readRDS(file = "output/NMF_comb_signature.rds")
show_sig_profile(Sig.Comb)
