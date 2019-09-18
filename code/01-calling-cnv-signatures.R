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
  pre_wes.prob <- derive(cp_wes, cores = ncores, nrep = 3)
)

system.time(
  pre_wes.count <- derive(cp_wes, type = "count", cores = ncores, nrep = 3)
)

save(pre_wes.prob, file = "output/pre_wes.prob.RData")
save(pre_wes.count, file = "output/pre_wes.count.RData")

sig.prob = sig_auto_extract(pre_wes.prob$nmf_matrix, result_prefix = "BayesNMF_Prob",
                            destdir = "output/signature", cores = 10)
sig.count = sig_auto_extract(pre_wes.count$nmf_matrix, result_prefix = "BayesNMF_Count",
                             destdir = "output/signature", cores = 10)
# BayesNMF will get more signatures and sparse results and
# it is hard to explain in this project

est_wes.prob = sig_estimate(pre_wes.prob$nmf_matrix,
                       range = 2:10, nrun = 50, cores = ncores, use_random = TRUE,
                       save_plots = TRUE, plot_basename = "output/plots/plot_estimate/wes_prob",
                       verbose = TRUE)
est_wes.count = sig_estimate(pre_wes.count$nmf_matrix,
                            range = 2:10, nrun = 50, cores = ncores, use_random = TRUE,
                            save_plots = TRUE, plot_basename = "output/plots/plot_estimate/wes_count",
                            verbose = TRUE)

save(est_wes.prob, file = "output/est_wes.prob.RData")
save(est_wes.count, file = "output/est_wes.count.RData")


rank_sury = rbind(est_wes.count$survey)
rank_sury = rbind(est_wes.prob$survey)

par(mar = c(5, 5, 2, 6))
plot(rank_sury$rank, rank_sury$cophenetic, type = "b", ann = FALSE)
mtext("Stability (cophenetic)", side = 2, line = 3)
mtext("Number of signature", side = 1, line = 3)
par(new = TRUE)
plot(x = rank_sury$rank, y = rank_sury$rss, type = "b", col = "red", ann = FALSE, axes = FALSE)
mtext("Error (RSS)", side = 4, line = 3, las = 0, col = "red")
axis(4)


# Use NMF instead of bayesian NMF
sig.prob2 = sig_extract(pre_wes.prob$nmf_matrix, n_sig = 5, nrun = 50, cores = ncores)
sig.count2 = sig_extract(pre_wes.count$nmf_matrix, n_sig = 5, nrun = 50, cores = ncores)

saveRDS(sig.prob2, file = "output/NMF_copynumber_signature.prob.rds")
saveRDS(sig.count2, file = "output/NMF_copynumber_signature.count.rds")


# Normalise by Row
show_sig_profile(sig.prob2, params = pre_wes.prob$parameters, y_expand = 1.5)
show_sig_profile(sig.count2, params = pre_wes.count$parameters, y_expand = 1.5)

# Normalise by Column
show_sig_profile(sig.prob2, params = pre_wes.prob$parameters, y_expand = 1.5, normalize = "column")
show_sig_profile(sig.count2, params = pre_wes.count$parameters, y_expand = 1.5, normalize = "column")




# cairo_pdf("wgs/plot_sig/cn_signature.pdf", height = 5)
# draw_sig_profile(sigs$nmfObj, params = h_data$parameters, y_expand = 2.5, sig_names = paste0("Sig", 1:4))
# dev.off()
# cowplot::ggsave2("wgs/plot_sig/cn_signature.pdf")
