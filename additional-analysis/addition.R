library(sigminer)
load("output/CNV.seqz.RData")
df.seqz <- readRDS("output/df.seqz.RDS")
tcga_prad <- readRDS("additional-analysis/prad.rds")
tcga_ov <- readRDS("additional-analysis/OV.rds")

# Compare CNA across cohorts ---------------------------------
# 包括 profile 汇总，signature 贡献
show_group_distribution(data = df.seqz %>%
                          dplyr::select(c("Study", "n_CNV")) %>%
                          na.omit, gvar = "Study", dvar = "n_CNV", order_by_fun = TRUE, ylab = "Number of CNA")

table(df.seqz$Study, df.seqz$sample_type)

show_group_distribution(data = df.seqz %>%
                          dplyr::select(c("sample_type", "CN-Sig2")) %>%
                          na.omit, gvar = "sample_type", dvar = "CN-Sig2", order_by_fun = TRUE, ylab = "Exposure of CN-Sig2")

# Chromothripsis score checking --------------------------------------------------

library(sigminer)

score_df = scoring(CNV.seqz, TD_size_cutoff = c(1e3, 1e5, 2e6))
score_df

plot_circos = function(object, data, col, text = "", decreasing = TRUE, median = FALSE, score = NULL,
                       color = circlize::colorRamp2(c(1, 2, 4), c("blue", "black", "red"))) {

  if (is.null(score)) {
    if (median) {
      samp = data[order(data[[col]],
                        decreasing = decreasing)]$sample[round(length(data$sample)/2)]
    } else {
      samp = data[order(data[[col]],
                        decreasing = decreasing)]$sample[1]
    }
  } else {
    samp = data$sample[data[[col]] == score][1]
  }
  print(samp)
  score = data[[col]][data$sample == samp]
  print(data[data$sample == samp])
  show_cn_circos(object, samples = samp,
                 col = color)
  text(0, 0, paste0(text, "\nscore = ", round(score, 3)))
}



pdf("additional-analysis/ChromothripisisScore-Check.pdf", width = 8, height = 8.5)
layout(matrix(1:4, nrow=2, ncol=2, byrow = TRUE))
plot_circos(CNV.seqz, score_df, col = "Chromoth_state", text = "Chromothripisis", decreasing = FALSE)
plot_circos(CNV.seqz, score_df, col = "Chromoth_state", text = "Chromothripisis", median = TRUE)
plot_circos(CNV.seqz, score_df, col = "Chromoth_state", text = "Chromothripisis", score = 86)
plot_circos(CNV.seqz, score_df, col = "Chromoth_state", text = "Chromothripisis")
dev.off()
layout(1)



# Validate WES cn profile with TCGA cn profile ----------------------------
tcga_prad$sample = substr(tcga_prad$sample_submitter_id, 1, 15)

# Set this per R session
options(sigminer.sex = "male", sigminer.copynumber.max = 20L)

# Generate CopyNumber object ----------------------------------------------

CNV.tcga_prad <- read_copynumber(tcga_prad,
                                 seg_cols = c("Chromosome", "Start", "End", "Copy_Number"),
                                 genome_build = "hg38",
                                 complement = FALSE,
                                 verbose = TRUE
)

CNV.tcga_prad_wes <- subset(CNV.seqz, startsWith(sample, "TCGA"))

common_samps <- CNV.tcga_prad@summary.per.sample$sample[CNV.tcga_prad@summary.per.sample$sample %in% CNV.tcga_prad_wes@summary.per.sample$sample]
CNV.tcga_prad <- subset(CNV.tcga_prad, sample %in% common_samps)
CNV.tcga_prad_wes <- subset(CNV.tcga_prad_wes, sample %in% common_samps)

save(CNV.tcga_prad, CNV.tcga_prad_wes, file = "additional-analysis/CNV_TCGA_PRAD.RData")

df1 = data.table::copy(CNV.tcga_prad@data)
df2 = data.table::copy(CNV.tcga_prad_wes@data)

# Calculate the summary that the segments can be roughly matched.
cos <- sapply(common_samps, function(id) {
  dt1 <- df1[sample == id]
  dt2 <- df2[sample == id]
  data.table::setkey(dt2, chromosome, start, end)
  y = data.table::foverlaps(dt1, dt2)
  y = na.omit(y)
  sigminer::cosine(y$segVal, y$i.segVal)
})

hist(cos, xlab = NA, main = "Cosine similarity of\nTCGA PRAD CN data from WES and SNP")

tail(sort(cos), 50)

library(ggplot2)

p11 <- show_cn_profile(CNV.tcga_prad, samples = "TCGA-EJ-A46D-01", show_title = TRUE)
p12 <- show_cn_profile(CNV.tcga_prad_wes, samples = "TCGA-EJ-A46D-01", ylim = c(0, 6))

p21 <- show_cn_profile(CNV.tcga_prad, samples = "TCGA-KK-A59X-01", ylim = c(0, 4), show_title = TRUE)
p22 <- show_cn_profile(CNV.tcga_prad_wes, samples = "TCGA-KK-A59X-01")

p31 <- show_cn_profile(CNV.tcga_prad, samples = "TCGA-EJ-7784-01", ylim = c(2, 7), show_title = TRUE)
p32 <- show_cn_profile(CNV.tcga_prad_wes, samples = "TCGA-EJ-7784-01")

library(patchwork)
pdf("additional-analysis/WES-SNP-CN-match.pdf", width = 15, height = 7)
old_par <- par(mar = c(2, 2.2, 0, 0), bg = NA)
wrap_elements(panel = ~hist(cos, xlab = NA, breaks = 30, main = "Cosine similarity of\nTCGA PRAD CN data from WES and SNP"), clip = FALSE) /
  ((p11 / p12) | (p21 / p22) | (p31 / p32))
par(old_par)
dev.off()

# 486 samples
library(NMF)
tally_W = sig_tally(CNV.tcga_prad, method = "W", cores = 6)
Sig.PRAD_SNP <- sig_auto_extract(tally_W$nmf_matrix, nrun = 5, cores = 5, K0 = 5)

load("output/Sig.CNV.seqz.W.RData")

get_sig_similarity(Sig.CNV.seqz.W, Sig.PRAD_SNP, normalize = "feature",  pattern_to_rm = "BoChr")

# An show case to other cancer types and fixed components -----------------

options(sigminer.sex = "female", sigminer.copynumber.max = 20L)

# Generate CopyNumber object ----------------------------------------------

CNV.tcga_ov <- read_copynumber(tcga_ov,
                               seg_cols = c("Chromosome", "Start", "End", "Copy_Number"),
                               genome_build = "hg38",
                               complement = FALSE,
                               verbose = TRUE
)

library(NMF)
tally_W_ov = sig_tally(CNV.tcga_ov, method = "W", cores = 6)
Sig.OV_SNP <- sig_auto_extract(tally_W$nmf_matrix, nrun = 10, cores = 5, K0 = 15)
