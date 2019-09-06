library(sigminer)
library(NMF)
library(data.table)


# Reading data ------------------------------------------------------------

plan_read = drake_plan(
  seg_wes = rbind(fread(file_in("data/CNV_from_TCGA_WES.tsv")),
                  fread(file_in("data/CNV_from_dbGap_except_TCGA_WES.tsv"))),
  cp_wes = read_copynumber(input = seg_wes, genome_build = "hg38", complement = FALSE, verbose = TRUE),
  seg_wgs = fread(file_in("data/CNV_from_TCGA_WGS.tsv")),
  cp_wgs = read_copynumber(input = seg_wgs, genome_build = "hg19", complement = FALSE, verbose = TRUE)
)

vis_plan(plan_read, font_size = 15)
make(plan_read, jobs = 2, parallelism = "future")

# Load data from cache
loadd(cp_wes, cp_wgs)

# Check distribution
draw_cn_distribution(cp_wes)
draw_cn_distribution(cp_wes, mode = "cd", fill = TRUE)

draw_cn_distribution(cp_wgs)
draw_cn_distribution(cp_wgs, mode = "cd", fill = TRUE)

boxplot(cp_wes@summary.per.sample$cna_burden)


# Prepare data and estimate rank -------------------------------------------

ncores = 16

#pre_wes = sig_prepare(cp_wes, cores = ncores)
plan_prepare = drake_plan(
  pre_wes = sig_prepare(cp_wes, cores = ncores),
  est_wes = sig_estimate(pre_wes$nmf_matrix,
                         range = 2:20, nrun = 100, cores = ncores, use_random = TRUE,
                         save_plots = TRUE, plot_basename = "output/plots/plot_estimate/wes",
                         verbose = TRUE),
  pre_wgs = sig_prepare(cp_wgs, reference_components = pre_wes$components, cores = ncores)
)

vis_plan(plan_prepare)
make(plan_prepare, lock_envir = FALSE)

loadd(pre_wes, pre_wgs)
#loadd(est_wes)

rank_sury = rbind(est_wes$survey)

par(mar = c(5, 5, 2, 6))
plot(rank_sury$rank, rank_sury$cophenetic, type = "b", ann = FALSE)
mtext("Stability (cophenetic)", side = 2, line = 3)
mtext("Number of signature", side = 1, line = 3)
par(new = TRUE)
plot(x = rank_sury$rank, y = rank_sury$rss, type = "b", col = "red", ann = FALSE, axes = FALSE)
mtext("Error (RSS)", side = 4, line = 3, las = 0, col = "red")
axis(4)
# select 5 or 6

# Choose best signature number
sigs5 = sig_extract(pre_wes$nmf_matrix, n_sig = 5, nrun = 100, cores = ncores)
sigs6 = sig_extract(pre_wes$nmf_matrix, n_sig = 6, nrun = 100, cores = ncores)

draw_cn_features(pre_wes$features)
params = draw_cn_components(pre_wes$features, pre_wes$components)
draw_sig_profile(sigs5$nmfObj, params = params$parameters, y_expand = 3)
draw_sig_profile(sigs6$nmfObj, params = params$parameters, y_expand = 3)

draw_cn_features(cn.pre$features)
draw_cn_components(cn.pre$features, cn.pre$components)

estimate = sig_estimate(cn.pre$nmf_matrix,
                        )
sigs = sig_extract(cn.pre$nmf_matrix, n_sig = 4, nrun = 100, cores = 4)
h_data = draw_cn_components(cn.pre$features, cn.pre$components)

cairo_pdf("wgs/plot_sig/cn_signature.pdf", height = 5)
draw_sig_profile(sigs$nmfObj, params = h_data$parameters, y_expand = 2.5, sig_names = paste0("Sig", 1:4))
dev.off()
cowplot::ggsave2("wgs/plot_sig/cn_signature.pdf")
