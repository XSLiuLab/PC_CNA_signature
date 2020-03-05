library(sigminer)
load("output/CNV.seqz.RData")
df = readRDS("output/df.seqz.RDS")

score_df = scoring(CNV.seqz, TD_size_cutoff = c(1e3, 1e5, 2e6))
score_df

tb = merge(df[, c("CNV_ID", "cnv_enrich_sig")], score_df, by.x = "CNV_ID", by.y = "sample")
boxplot(MACN ~ cnv_enrich_sig, data = tb)

# boxplot(sTD ~ cnv_enrich_sig, data = tb)
#
# boxplot(TD ~ cnv_enrich_sig, data = tb)
# boxplot(lTD ~ cnv_enrich_sig, data = tb)
#
# boxplot(sTDP ~ cnv_enrich_sig, data = tb)
# boxplot(lTDP ~ cnv_enrich_sig, data = tb)
#
# tb_cnv = CNV.seqz@data[segVal > 2]
# tb_cnv$ss = log10(tb_cnv$end - tb_cnv$start) + 1
# hist(tb_cnv$ss, breaks = 100)
# plot(density(tb_cnv$ss), xlim = c(3, 8))
# abline(v=5, lty = 2)
# abline(v=6, lty = 2)
# abline(v=7, lty = 2)
#
# #Sig1 enriched
# plot(density(tb_cnv$ss[tb_cnv$sample == "5115056-SRR8311885" & tb_cnv$segVal > 8]), xlim = c(3, 8),
#      main = "5115056-SRR8311885")
# hist(tb_cnv$ss[tb_cnv$sample == "5115056-SRR8311885" & tb_cnv$segVal > 8], breaks = 100)
#
#
# #Sig2 enriched
# plot(density(tb_cnv$ss[tb_cnv$sample == "5115615-SRR8311749" & tb_cnv$segVal <= 4 & tb_cnv$segVal > 2]), xlim = c(3, 8),
#      main = "5115615-SRR8311749")
#
# hist(tb_cnv$ss[tb_cnv$sample == "5115615-SRR8311749" & tb_cnv$segVal <= 4 & tb_cnv$segVal > 2],
#      breaks = 100)
#
# #Sig3 enriched
# plot(density(tb_cnv$ss[tb_cnv$sample == "TCGA-KK-A7B4-01"]), xlim = c(3, 8),
#      main = "TCGA-KK-A7B4-01")
#
# hist(tb_cnv$ss[tb_cnv$sample == "TCGA-KK-A7B4-01"],
#      breaks = 100)



## check
pdf("score_density_checking.pdf", width = 7, height = 6)
# layout(matrix(1:9, nrow=3, ncol=3, byrow = TRUE))
layout(matrix(c(1:5, 0, 6, 0, 0), nrow=3, ncol=3, byrow = TRUE))

# plot(density(score_df$TDP, na.rm = T), main = "TDP score")
plot(density(score_df$TD, na.rm = T), main = "TD score")
# plot(density(score_df$sTD, na.rm = T), main = "short TD score")
# plot(density(score_df$lTD, na.rm = T), main = "long TD score")
plot(density(score_df$cnaBurden, na.rm = T), main = "CNA burden score")
# plot(density(score_df$cnaLoad, na.rm = T), main = "CNA load score")
plot(density(score_df$MACN, na.rm = T), main = "MACN score")
# plot(density(score_df$weightedMACN, na.rm = T), main = "Weighted MACN score")
plot(density(score_df$Ploidy, na.rm = T), main = " ploidy score")
plot(density(log2(score_df$Chromothripisis), na.rm = T), "log2 Chromothripisis score")

plot(score_df$TDP, score_df$TD,
     xlab = "TDP score (method from PNAS)",
     ylab = "TD score (our definition)",
     pch = 2, cex = 0.5,
     main = paste("Correlation coefficient r =",
                  round(cor(score_df$TDP, score_df$TD,
                            use = "pairwise.complete.obs",
                            method = "spearman"), 3)))
dev.off()
layout(1)


plot_circos = function(object, data, col, text = "", decreasing = TRUE,
                       color = circlize::colorRamp2(c(1, 2, 4), c("blue", "white", "red"))) {

  samp = data[order(data[[col]],
                      decreasing = decreasing)]$sample[1]
  print(samp)
  score = data[[col]][data$sample == samp]
  print(data[data$sample == samp])
  show_cn_circos(object, samples = samp,
                 col = color)
  text(0, 0, paste0(text, "\nscore = ", round(score, 3)))
}


pdf("score_circos_heatmap.pdf", width = 10, height = 8.5)
layout(matrix(1:12, nrow=3, ncol=4, byrow = TRUE))
plot_circos(CNV.seqz, score_df, col = "cnaBurden", text = "CNA burden",
            color = circlize::colorRamp2(c(1, 2, 6), c("blue", "white", "red")))
plot_circos(CNV.seqz, score_df, col = "cnaBurden", text = "CNA burden", decreasing = FALSE)
plot_circos(CNV.seqz, score_df, col = "cnaLoad", text = "CNA number")
plot_circos(CNV.seqz, score_df, col = "cnaLoad", text = "CNA number", decreasing = FALSE)

# cor(score_df$MACN, score_df$weightedMACN, use = "pairwise.complete.obs")
plot_circos(CNV.seqz, score_df, col = "MACN", text = "MACN")
plot_circos(CNV.seqz, score_df, col = "MACN", text = "MACN", decreasing = FALSE)
# plot_circos(CNV.seqz, score_df, col = "weightedMACN", text = "Weighted MACN")
# plot_circos(CNV.seqz, score_df, col = "weightedMACN", text = "Weighted MACN", decreasing = FALSE)

plot_circos(CNV.seqz, score_df, col = "Ploidy", text = "Ploidy",
            color = circlize::colorRamp2(c(1, 2, 6), c("blue", "white", "red")))
plot_circos(CNV.seqz, score_df, col = "Ploidy", text = "Ploidy", decreasing = FALSE)

plot_circos(CNV.seqz, score_df, col = "Chromothripisis", text = "Chromothripisis")
plot_circos(CNV.seqz, score_df, col = "Chromothripisis", text = "Chromothripisis", decreasing = FALSE)

plot_circos(CNV.seqz, score_df, col = "TD", text = "TD")
plot_circos(CNV.seqz, score_df, col = "TD", text = "TD", decreasing = FALSE)
# plot_circos(CNV.seqz, score_df, col = "sTD", text = "short TD")
# plot_circos(CNV.seqz, score_df, col = "sTD", text = "short TD", decreasing = FALSE)
# plot_circos(CNV.seqz, score_df, col = "lTD", text = "long TD")
# plot_circos(CNV.seqz, score_df, col = "lTD", text = "long TD", decreasing = FALSE)
dev.off()
layout(1)
