library(sigminer)
library(NMF)
library(data.table)

# Detecting copy number signatures ----------------------------------------

segtab <- fread("wgs/TCGA_WGS_facets_CNV.txt")
cn <- read_copynumber(input = segtab, genome_build = "hg19", complement = FALSE, verbose = TRUE)

# Check distribution
draw_cn_distribution(cn)
draw_cn_distribution(cn, mode = "cd")
draw_cn_distribution(cn, mode = "cd", fill = TRUE)

cn.pre <- sig_prepare(cn, cores = 4, seed = 123456)

draw_cn_features(cn.pre$features)
draw_cn_components(cn.pre$features, cn.pre$components)

estimate <- sig_estimate(cn.pre$nmf_matrix,
  range = 2:10, nrun = 100, cores = 4, use_random = TRUE,
  save_plots = TRUE, plot_basename = "wgs/plot_estimate/tcga_wgs",
  verbose = TRUE
)
sigs <- sig_extract(cn.pre$nmf_matrix, n_sig = 4, nrun = 100, cores = 4)
h_data <- draw_cn_components(cn.pre$features, cn.pre$components)

cairo_pdf("wgs/plot_sig/cn_signature.pdf", height = 5)
draw_sig_profile(sigs$nmfObj, params = h_data$parameters, y_expand = 2.5, sig_names = paste0("Sig", 1:4))
dev.off()
cowplot::ggsave2("wgs/plot_sig/cn_signature.pdf")
