load("output/Sig.PRAD_TCGA_plus_dbGap_Maf.RData")
load("output/Sig.CNV.seqz.W.RData")

NMF::consensusmap(Sig.CNV.seqz.W$Raw$nmf_obj,
                  tracks =c("consensus:", "silhouette:"))

png(file = "consensus_map_snv.png", width = 5, height = 4, units = "in", res = 300)
NMF::consensusmap(Sig.SNV$Raw$nmf_obj,
                  tracks =c("consensus:", "silhouette:"),
                  labCol = NA, labRow = NA)
dev.off()
