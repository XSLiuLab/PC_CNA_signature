library(sigminer)
#library(drake)
library(NMF)
library(data.table)
library(dtplyr)


# Reading data ------------------------------------------------------------

CNV = data.table::fread("data/CNV_from_dbGAP_PLUS_TCGA_WES_CVAL150.tsv")
CNV_X = CNV[Chromosome == 23]
CNV = CNV[Chromosome != 23]
# Double the copy number of X due to all samples are male
CNV_X$modal_cn = 2 * CNV_X$modal_cn
CNV = rbind(CNV, CNV_X)

# Check number of segment for samples
CNV %>%
  lazy_dt() %>%
  group_by(sample) %>%
  summarize(n_chr = length(unique(Chromosome))) %>%
  arrange(n_chr) %>%
  data.table::as.data.table()
# find 141-10 only have segments in chr1
# remove it
CNV = CNV[!sample %in% "141-10"]

table(CNV$Chromosome)
CNV = read_copynumber(CNV, genome_build = "hg38",
                      complement = FALSE, verbose = TRUE)
save(CNV, file="data/PRAD_CNV.RData")

CNV@summary.per.sample$sample[startsWith(CNV@summary.per.sample$sample, prefix = "447")]

hist(CNV@data$segVal, breaks = 20)

# Check distribution
show_cn_distribution(CNV)
show_cn_distribution(CNV, mode = "cd", fill = TRUE)

boxplot(CNV@summary.per.sample$cna_burden)


# Prepare data and estimate rank -------------------------------------------
ncores = 12

system.time(
  CNV.prob <- sig_derive(CNV, cores = ncores, nrep = 3)
)

system.time(
  CNV.count <- sig_derive(CNV, type = "count", cores = ncores, nrep = 3)
)

save(CNV.prob, file = "output/CNV.prob.RData")
save(CNV.count, file = "output/CNV.count.RData")

show_cn_features(CNV.prob$features)
show_cn_features(CNV.prob$features, log_segsize = F)
show_cn_components(CNV.prob$parameters)
show_cn_components(CNV.prob$parameters, show_weights = FALSE, log_segsize = F)

# Auto-extract
Sig.Bayesian.prob = sig_auto_extract(CNV.prob$nmf_matrix, result_prefix = "BayesNMF_Prob", nrun = 100,
                                     destdir = "output/signature", cores = 16)
Sig.Bayesian.count = sig_auto_extract(CNV.count$nmf_matrix, result_prefix = "BayesNMF_Count", nrun = 100,
                                     destdir = "output/signature", cores = 16)

# BayesNMF will get more signatures and sparse results and
# it is hard to explain in this project

CNV_EST.prob = sig_estimate(CNV.prob$nmf_matrix,
                       range = 2:10, nrun = 50, cores = ncores, use_random = TRUE,
                       save_plots = FALSE,
                       verbose = TRUE)
CNV_EST.count = sig_estimate(CNV.count$nmf_matrix,
                            range = 2:10, nrun = 50, cores = ncores, use_random = TRUE,
                            save_plots = FALSE ,pConstant = 0.001,
                            verbose = TRUE)

save(CNV_EST.prob, file = "output/CNV_EST.prob.RData")
save(CNV_EST.count, file = "output/CNV_EST.count.RData")

show_rank_survey(CNV_EST.prob)
show_rank_survey(CNV_EST.count)

# Use NMF instead of bayesian NMF
Sig.CNV.prob = sig_extract(CNV.prob$nmf_matrix, n_sig = 6, nrun = 100, cores = ncores)
Sig.CNV.count = sig_extract(CNV.count$nmf_matrix, n_sig = 6, nrun = 100, cores = ncores, pConstant = 0.001)

saveRDS(Sig.CNV.prob, file = "output/NMF_copynumber_signature.prob.rds")
saveRDS(Sig.CNV.count, file = "output/NMF_copynumber_signature.count.rds")

# Take a look
# Normalise by Row
show_sig_profile(Sig.CNV.prob, params = CNV.prob$parameters, y_expand = 1.5)
show_sig_profile(Sig.CNV.count, params = CNV.count$parameters, y_expand = 1.5)

# Normalise by Column
show_sig_profile(Sig.CNV.prob, params = CNV.prob$parameters, y_expand = 1.5, normalize = "column")
show_sig_profile(Sig.CNV.count, params = CNV.count$parameters, y_expand = 1.5, normalize = "column")


# cairo_pdf("wgs/plot_sig/cn_signature.pdf", height = 5)
# draw_sig_profile(sigs$nmfObj, params = h_data$parameters, y_expand = 2.5, sig_names = paste0("Sig", 1:4))
# dev.off()
# cowplot::ggsave2("wgs/plot_sig/cn_signature.pdf")
