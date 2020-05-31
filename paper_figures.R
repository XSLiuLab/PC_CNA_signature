## Generate paper figures
## Most of figures are directly generated from RMarkdown files

load("output/Sig.PRAD_TCGA_plus_dbGap_Maf.RData")
load("output/Sig.CNV.seqz.W.RData")


# Main --------------------------------------------------------------------

## Figure 2
p = show_sig_profile(Sig.SNV, mode = "SBS", style = "cosmic", x_label_angle = 90, x_label_vjust = 0.5)
p = add_labels(p, x = 0.92, y = 0.3, y_end = 0.85, n_label = 3,
           labels = rev(c("HRD or unknown", "dMMR", "Aging")), hjust = 1)

ggsave("figures/Figure2A.pdf", width = 10, height = 5)


p = show_sig_profile(Sig.CNV.seqz.W, style = "cosmic", method = "W", normalize = "feature")
ggsave("figures/Figure2B.pdf", width = 12, height = 7)


# Supp figures ------------------------------------------------------------


NMF::consensusmap(Sig.CNV.seqz.W$Raw$nmf_obj,
                  tracks =c("consensus:", "silhouette:"))

png(file = "figures/Figure_S4_C.png", width = 5, height = 4, units = "in", res = 300)
NMF::consensusmap(Sig.CNV.seqz.W$Raw$nmf_obj,
                  tracks =c("consensus:", "silhouette:"),
                  labCol = NA, labRow = NA)
dev.off()

png(file = "figures/Figure_S4_D.png", width = 5, height = 4, units = "in", res = 300)
NMF::consensusmap(Sig.SNV$Raw$nmf_obj,
                  tracks =c("consensus:", "silhouette:"),
                  labCol = NA, labRow = NA)
dev.off()
