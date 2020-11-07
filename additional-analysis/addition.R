library(sigminer)
library(ggplot2)
load("output/CNV.seqz.RData")
df.seqz <- readRDS("output/df.seqz.RDS")
tcga_prad <- readRDS("additional-analysis/prad.rds")
tcga_ov <- readRDS("additional-analysis/OV.rds")

# Show CNA frequency plot for different studies ---------------------------
df.seqz %>%
  dplyr::select(Study, CNV_ID) %>%
  dplyr::group_split(Study) %>%
  purrr::map(.f = function(df) {
    study = df$Study[1]
    samps = df$CNV_ID
    obj <- subset(CNV.seqz, sample %in% samps)
    p <- show_cn_group_profile(obj) + ggplot2::labs(title = study)
  }) -> p_list

p <- cowplot::plot_grid(plotlist = p_list[c(6,1:5)], align = "hv", ncol = 2)
ggplot2::ggsave(filename = "additional-analysis/study_group_profile.pdf",
                plot = p,
                width = 10,
                height = 8)

# Compare CNA across cohorts ---------------------------------
# 包括 profile 汇总，signature 贡献
p1 <- show_group_distribution(data = df.seqz %>%
                          dplyr::select(c("Study", "n_CNV")) %>%
                          na.omit, gvar = "Study", dvar = "n_CNV", order_by_fun = TRUE, ylab = "Number of CNA")
ggsave("additional-analysis/study-cna-dist.pdf", plot = p1, width = 6, height = 3)

table(df.seqz$Study, df.seqz$sample_type)

p2 <- show_group_distribution(data = df.seqz %>%
                          dplyr::select(c("sample_type", "CN-Sig2")) %>%
                          na.omit, gvar = "sample_type", dvar = "CN-Sig2", order_by_fun = TRUE, ylab = "Exposure of CN-Sig2")
# show it or not?

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
                               samp_col = "sample_submitter_id",
                               genome_build = "hg38",
                               complement = FALSE,
                               verbose = TRUE
)

library(NMF)

# Pick up 50 samples
# generate subset 10, 25, and 50
# and see how method "W" and "M" would work for the components
set.seed(1234)
test_samples = sample(CNV.tcga_ov@summary.per.sample$sample, 50)
test_samples25 = sample(test_samples, 25)
test_samples10 = sample(test_samples25, 10)

test_list <- list(test_samples, test_samples25, test_samples10)

tally_W_list <- list()
tally_M_list <- list()
for (i in seq_along(test_list)) {
  CNV.test <- subset(CNV.tcga_ov, sample %in% test_list[[i]])
  tally_W_list[[i]] = sig_tally(CNV.test, method = "W", cores = 6)
  tally_M_list[[i]] = sig_tally(CNV.test, method = "M", cores = 6)
}


tally_M_list[[1]]$nmf_matrix %>% dim()
tally_M_list[[2]]$nmf_matrix %>% dim()
tally_M_list[[3]]$nmf_matrix %>% dim()

est_W <- sig_estimate(tally_W_list[[1]]$nmf_matrix, range = 2:10, nrun = 30, pConstant = 1e-9, verbose = TRUE)
show_sig_number_survey(est_W$survey)

est_M <- sig_estimate(tally_M_list[[1]]$nmf_matrix, range = 2:10, nrun = 30, pConstant = 1e-9, verbose = TRUE)
show_sig_number_survey(est_M$survey)

# Just take 4 signatures
sigs_W_list <- list()
sigs_M_list <- list()
for (i in seq_along(tally_W_list)) {
  sigs_W_list[[i]] = sig_extract(tally_W_list[[i]]$nmf_matrix, n_sig = 4, nrun = 30, cores = 6, pConstant = 1e-9)
  sigs_M_list[[i]] = sig_extract(tally_M_list[[i]]$nmf_matrix, n_sig = 4, nrun = 30, cores = 6, pConstant = 1e-9)
}

p11 <- show_sig_profile(sigs_W_list[[1]], mode = "copynumber", method = "W", normalize = "feature", style = "cosmic")
p12 <- show_sig_profile(sigs_W_list[[2]], mode = "copynumber", method = "W", normalize = "feature", style = "cosmic")
p13 <- show_sig_profile(sigs_W_list[[3]], mode = "copynumber", method = "W", normalize = "feature", style = "cosmic")

get_sig_similarity(sigs_W_list[[2]], sigs_W_list[[1]])
get_sig_similarity(sigs_W_list[[3]], sigs_W_list[[1]])

p21 <- show_sig_profile(sigs_M_list[[1]], mode = "copynumber", method = "M", params = tally_M_list[[1]]$parameters,
                 normalize = "feature", style = "cosmic", paint_axis_text = F, y_expand = 1.5)
p22 <- show_sig_profile(sigs_M_list[[2]], mode = "copynumber", method = "M", params = tally_M_list[[2]]$parameters,
                 normalize = "feature", style = "cosmic", paint_axis_text = F, y_expand = 1.5)
p23 <- show_sig_profile(sigs_M_list[[3]], mode = "copynumber", method = "M", params = tally_M_list[[3]]$parameters,
                 normalize = "feature", style = "cosmic", paint_axis_text = F, y_expand = 1.5)

library(patchwork)
p1 <- (p11 / p12 / p13)
p2 <- (p21 / p22 / p23)

ggsave("additional-analysis/ovary_sigprofile_W.pdf", plot = p1, width = 14, height = 12)
ggsave("additional-analysis/ovary_sigprofile_M.pdf", plot = p2, width = 6, height = 12, device = cairo_pdf)

sim1 <- get_sig_similarity(sigs_W_list[[2]], sigs_W_list[[1]])
sim2 <- get_sig_similarity(sigs_W_list[[3]], sigs_W_list[[1]])

pheatmap::pheatmap(sim1$similarity, display_numbers = TRUE,
                   width = 3, height = 3,
                   filename = "additional-analysis/25-vs-50-similarity.pdf")
pheatmap::pheatmap(sim2$similarity, display_numbers = TRUE,
                   width = 3, height = 3,
                   filename = "additional-analysis/10-vs-50-similarity.pdf")
