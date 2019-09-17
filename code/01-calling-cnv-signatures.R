library(sigminer)
library(drake)
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

hist(cp_wes@data$segVal, breaks = 20)
hist(cp_wgs@data$segVal, breaks = 20)
sort(cp_wes@data$segVal, decreasing = TRUE) %>% head(200)
sort(cp_wgs@data$segVal, decreasing = TRUE) %>% head(200)

# Check distribution
draw_cn_distribution(cp_wes)
draw_cn_distribution(cp_wes, mode = "cd", fill = TRUE)

draw_cn_distribution(cp_wgs)
draw_cn_distribution(cp_wgs, mode = "cd", fill = TRUE)

boxplot(cp_wes@summary.per.sample$cna_burden)


# Prepare data and estimate rank -------------------------------------------
ncores = 8

system.time(
  pre_wes <- derive(cp_wes, cores = ncores, nrep = 3)
)

system.time(
  pre_wes2 <- derive(cp_wes, cores = ncores, nrep = 3)
)


#pre_wes = sig_prepare(cp_wes, cores = ncores)
plan_prepare = drake_plan(
  pre_wes = derive(cp_wes, cores = ncores, nrep = 3),
  pre_wgs = derive(cp_wgs, reference_components = pre_wes$components, cores = ncores)
)

system.time(
  pre_wgs <- sig_prepare(cp_wgs, cores = ncores, nrep = 1)
)



save(pre_wes, file = "output/pre_wes.RData")

vis_plan(plan_prepare)
make(plan_prepare, lock_envir = FALSE)

loadd(pre_wes, pre_wgs)

# est_wes = sig_estimate(pre_wes$nmf_matrix,
#                        range = 2:15, nrun = 100, cores = ncores, use_random = TRUE,
#                        save_plots = TRUE, plot_basename = "output/plots/plot_estimate/wes",
#                        verbose = TRUE),
# loadd(est_wes)

res = sig_auto_extract(pre_wes$nmf_matrix, result_prefix = "BayesNMF_WES", destdir = "output/BayesNMF/",
                 nrun = 200, K0 = 20, niter = 2e6, cores = ncores, skip = TRUE)

rank_sury = rbind(est_wes$survey)

par(mar = c(5, 5, 2, 6))
plot(rank_sury$rank, rank_sury$cophenetic, type = "b", ann = FALSE)
mtext("Stability (cophenetic)", side = 2, line = 3)
mtext("Number of signature", side = 1, line = 3)
par(new = TRUE)
plot(x = rank_sury$rank, y = rank_sury$rss, type = "b", col = "red", ann = FALSE, axes = FALSE)
mtext("Error (RSS)", side = 4, line = 3, las = 0, col = "red")
axis(4)
# select 6 or 8

# Choose best signature number
sigs = sig_extract(pre_wes$nmf_matrix, n_sig = 8, nrun = 100, cores = ncores)

draw_cn_features(pre_wes$features)
draw_cn_components(pre_wes$features, pre_wgs$components)

params = draw_cn_components(pre_wes$features, pre_wes$components)
draw_sig_profile(sigs$nmfObj, params = params$parameters, y_expand = 3)



estimate = sig_estimate(cn.pre$nmf_matrix,
                        )
sigs = sig_extract(cn.pre$nmf_matrix, n_sig = 4, nrun = 100, cores = 4)
h_data = draw_cn_components(cn.pre$features, cn.pre$components)

cairo_pdf("wgs/plot_sig/cn_signature.pdf", height = 5)
draw_sig_profile(sigs$nmfObj, params = h_data$parameters, y_expand = 2.5, sig_names = paste0("Sig", 1:4))
dev.off()
cowplot::ggsave2("wgs/plot_sig/cn_signature.pdf")
