library(sigminer)
library(tidyverse)
library(NMF)


# Read data ---------------------------------------------------------------


CNV_seqz <- data.table::fread("data/CNV_from_sequenza.tsv")

CNV_seqz %>%
  group_by(sample) %>%
  summarize(n_chr = length(unique(Chromosome))) %>%
  arrange(n_chr) %>%
  data.table::as.data.table()

# remove WCMC160-SRR3146971
CNV_seqz <- CNV_seqz[!sample %in% "WCMC160-SRR3146971"]


CNV_facets <- data.table::fread("data/CNV_from_dbGAP_PLUS_TCGA_WES_CVAL150.tsv")
CNV_facets %>%
  group_by(sample) %>%
  summarize(n_chr = length(unique(Chromosome))) %>%
  arrange(n_chr) %>%
  data.table::as.data.table()

CNV_facets <- CNV_facets[!sample %in% "141-10"]

CNV_seqz <- read_copynumber(CNV_seqz,
  genome_build = "hg38",
  complement = FALSE, verbose = TRUE
)

CNV_facets <- read_copynumber(CNV_facets,
  genome_build = "hg38",
  complement = FALSE, verbose = TRUE
)

save(CNV_seqz, CNV_facets, file = "data/PRAD_CNV_seqz_and_facets.RData")


# Prepare -----------------------------------------------------------------
load(file = "data/PRAD_CNV_seqz_and_facets.RData")

ncores <- 10

CNV.seqz.pre <- sig_derive(CNV_seqz, method = "W", feature_setting = CN.features, cores = ncores)
CNV.facets.pre <- sig_derive(CNV_facets, method = "W", feature_setting = CN.features, cores = ncores)

show_cn_features(CNV.seqz.pre$features, method = "W")
show_cn_features(CNV.facets.pre$features, method = "W")

show_cn_components(CNV.seqz.pre$parameters, method = "W")
show_cn_components(CNV.facets.pre$parameters, method = "W")

# seqz.sigs = sig_auto_extract(CNV.seqz.pre$nmf_matrix, result_prefix = "BayesNMF_seqz_wang", nrun = 100,
#                              destdir = "output/signature", cores = ncores, skip = TRUE)
# facets.sigs = sig_auto_extract(CNV.facets.pre$nmf_matrix, result_prefix = "BayesNMF_facets_wang", nrun = 100,
#                              destdir = "output/signature", cores = ncores)

seqz.est <- sig_estimate(CNV.seqz.pre$nmf_matrix,
  range = 2:10, nrun = 50, cores = ncores, use_random = FALSE,
  save_plots = FALSE,
  verbose = TRUE
)

show_rank_survey(seqz.est)

seqz.est.50 <- sig_estimate(CNV.seqz.pre$nmf_matrix[, 1:50],
  range = 2:10, nrun = 50, cores = ncores, use_random = FALSE,
  save_plots = FALSE,
  verbose = TRUE
)

show_rank_survey(seqz.est.50)

facets.est <- sig_estimate(CNV.facets.pre$nmf_matrix,
  range = 2:10, nrun = 50, cores = ncores, use_random = TRUE,
  save_plots = FALSE, pConstant = 0.001,
  verbose = TRUE
)

show_rank_survey(facets.est)

sigs.seqz <- sig_extract(CNV.seqz.pre$nmf_matrix, n_sig = 7, cores = 10)
sigs.seqz_5 <- sig_extract(CNV.seqz.pre$nmf_matrix, n_sig = 5, cores = 10)
sigs.seqz_50_5 <- sig_extract(CNV.seqz.pre$nmf_matrix[, 1:50], n_sig = 5, cores = 10)
sigs.seqz_50_7 <- sig_extract(CNV.seqz.pre$nmf_matrix[, 1:50], n_sig = 7, cores = 10)


save(CNV.seqz.pre, CNV.facets.pre, seqz.est, seqz.est.50, facets.est,
  sigs.seqz, sigs.seqz_50_5, sigs.seqz_50_7,
  file = "test.RData"
)

cc <- sig_auto_extract(CNV.seqz.pre$nmf_matrix, K0 = 7, cores = 10)
show_sig_profile(cc, method = "W", normalize = "column")

show_sig_profile(sigs.seqz, method = "W", normalize = "raw")
show_sig_profile(sigs.seqz, method = "W", normalize = "column")
show_sig_profile(sigs.seqz, method = "W", normalize = "feature")

show_sig_profile(sigs.seqz_5, method = "W", normalize = "row", x_label_angle = 90)
show_sig_profile(sigs.seqz_5, method = "W", normalize = "column", x_label_angle = 90)
show_sig_profile(sigs.seqz_5, method = "W", normalize = "feature")

sigs.seqz_5show_sig_exposure(sigs.seqz, rm_space = T)
show_sig_exposure(sigs.seqz_50_5, rm_space = T)

show_sig_profile(sigs.seqz_50_5, method = "W")
show_sig_profile(sigs.seqz_50_7, method = "W", normalize = "column")

get_groups(sigs.seqz_50_5, method = "consensus")$enrich_sig %>% table()
get_groups(sigs.seqz_50_7, method = "consensus")$enrich_sig %>% table()


CNV.seqz.pre$nmf_matrix[tt, 24] %>% sum()
CNV.seqz.pre$nmf_matrix[, 24] %>% sum()

sigs.facets <- sig_extract(CNV.seqz.pre$nmf_matrix, n_sig = 6, cores = 10)
show_sig_profile(sigs.facets, method = "W", x_label_angle = 90, normalize = "column")
show_sig_profile(sigs.facets, method = "W", x_label_angle = 90, normalize = "row")


# FACETS 6 signatures -----------------------------------------------------

load("test.RData")
load(file = "output/CNV.prob.RData")

Sig.M.facets <- readRDS(file = "output/NMF_copynumber_signature.prob.rds")
Sig.W.facets <- sig_extract(CNV.seqz.pre$nmf_matrix, n_sig = 6, cores = 10)

# Normalise by Row
show_sig_profile(Sig.M.facets, params = CNV.prob$parameters, y_expand = 1.5)
show_sig_profile(Sig.W.facets, method = "W", normalize = "row", x_label_angle = 90)

# Normalise by Column
show_sig_profile(Sig.M.facets, params = CNV.prob$parameters, y_expand = 1.5, normalize = "column")
show_sig_profile(Sig.W.facets, method = "W", normalize = "column", x_label_angle = 90)

# Normalize by feature
show_sig_profile(Sig.W.facets, method = "W", normalize = "feature", x_label_angle = 90)


show_sig_exposure(Sig.M.facets, rm_space = TRUE)
show_sig_exposure(Sig.W.facets, rm_space = TRUE)
