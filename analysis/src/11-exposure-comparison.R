library(sigminer)

load("output/Sig.CNV.seqz.W.RData")

load("output/Sig.CNV.mskcc.RData")
load("output/CNV.mskcc.tally.W.RData")

# Switch Sig1 and Sig2 to make the signature consistent
# See https://shixiangwang.github.io/prad_signature/#method-comparison-between-macintyre-et-al-and-our-study
Sig.CNV.mskcc = sig_modify_names(Sig.CNV.mskcc, new_names = paste0("Sig", c(2,1,3:5)))

expo_mskcc_actual = sig_fit(t(CNV.mskcc.tally.W$nmf_matrix), sig = Sig.CNV.mskcc,
                              mode = "copynumber",
                              return_class = "data.table")

colnames(expo_mskcc_actual) = c("sample", paste0("MSKCC-Sig", c(2, 1, 3:5)))

expo_mskcc_estimate = sig_fit(t(CNV.mskcc.tally.W$nmf_matrix), sig = Sig.CNV.seqz.W,
                              mode = "copynumber",
                              return_class = "data.table")

expo_df <- merge(expo_mskcc_actual, expo_mskcc_estimate, by = "sample", all = TRUE)

pdf(file = "figures/exposure_comparison_with_sig_fit.pdf", width = 6, height = 5.5)
pheatmap::pheatmap(cor(expo_df[, -1][, c(2,1,3:10)]), display_numbers = TRUE)
dev.off()

