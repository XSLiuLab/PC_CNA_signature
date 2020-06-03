library(sigminer)


# Bootstrap Analysis ------------------------------------------------------

load("output/CNV.seqz.tally.W.RData")
load("output/Sig.CNV.seqz.W.RData")

scale_labels <- ggplot2::scale_x_discrete(labels = paste("CN-Sig", 1:5))

bt_w <- sig_fit_bootstrap_batch(catalogue_matrix = CNV.seqz.tally.W$nmf_matrix %>% t(),
                                sig = Sig.CNV.seqz.W,
                                methods = "QP",
                                n = 1000,
                                p_val_thresholds = 0.01,
                                type = "relative",
                                mode = "copynumber",
                                job_id = "CNSig_W",
                                result_dir = "data/bootstrap",
                                use_parallel = TRUE)

save(bt_w, file = "data/bootstrap_w.RData")


load("data/bootstrap_w.RData")

p1 <- show_sig_bootstrap_stability(bt_w, add.params = list(alpha = 0.1), ylab = "Signature instability (MRSE)") +
  ggplot2::theme(legend.position = "none") +
  scale_labels
p2 <- show_sig_bootstrap_exposure(bt_w, add.params = list(alpha = 0.1), highlight = "red", sample = "5115615-SRR8311749") +
  ggplot2::theme(legend.position = "none") +
  scale_labels
p3 <- show_sig_bootstrap_exposure(bt_w, add.params = list(alpha = 0.1), highlight = "red", sample = "TCGA-KK-A7B4-01") +
  ggplot2::theme(legend.position = "none") +
  scale_labels
#show_sig_bootstrap_error(bt_w, add.params = list(alpha = 0.1), highlight = "red") + ggplot2::theme(legend.position = "none")


# Sampling analysis -------------------------------------------------------
# Samples a fraction of samples
# and calculate the prob a signature can be detected with
# cutoff at 1% relative exposure (p < 0.05) in at least 10 samples.

p_dt <- data.table::copy(bt_w$p_val)
p_dt <- p_dt[, .(isExist = p_value < 0.05), by = .(sample, sig)]

#fractions <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
fractions <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5)
samps <- unique(p_dt$sample)
n_samps <- length(samps)
res <- purrr::map_df(seq_along(fractions), function(i) {
  message("=> Starting for fraction ", fractions[i])
  n <- round(fractions[i] * n_samps)
  res <- purrr::map_df(1:1000,
                       .f = function(x) {
                         message("==> Doing run #", x)
                         subset_samps <- sample(samps, n, replace = FALSE)
                         df <- p_dt[sample %in% subset_samps][, .(Exist = sum(isExist) > 10), by = .(sig)]
                         dplyr::bind_cols(
                           dplyr::tibble(Run = x, .rows = nrow(df)),
                           df
                         )
                       }
  )
  dplyr::bind_cols(
    dplyr::tibble(Fraction = fractions[i], .rows = nrow(res)),
    res
  )
})

mp = paste("CN-Sig", 1:5)
names(mp) = paste0("Sig", 1:5)

res_summary <- res %>%
  dplyr::group_by(sig, Fraction) %>%
  dplyr::summarise(
    Prob = sum(Exist) / 1000
  ) %>%
  dplyr::ungroup() %>%
  dplyr::rename(
    Signature = sig
  ) %>%
  dplyr::mutate(
    Signature = mp[Signature] %>% as.character()
  )

library(ggplot2)
p4 <- ggplot(res_summary, aes(x = Fraction, y = Prob, color = Signature)) +
  geom_point() +
  geom_line() + cowplot::theme_cowplot() +
  geom_vline(xintercept = 0.06, linetype = 2) +
  labs(x = "Sampling fraction", y = "Detection probability") +
  theme(legend.position = "top") +
  scale_x_continuous(breaks = c(0, 0.06, 0.1, 0.2, 0.3, 0.4, 0.50)) +
  ggplot2::guides(color=ggplot2::guide_legend(nrow=2,byrow=TRUE))

937 * 0.06

# Fitting analysis --------------------------------------------------------

optim_w <- sig_fit(catalogue_matrix = CNV.seqz.tally.W$nmf_matrix %>% t(),
                   sig = Sig.CNV.seqz.W,
                   method = "QP",
                   type = "relative",
                   mode = "copynumber")

p5 <- show_sig_fit(optim_w, add.params = list(alpha = 0.1))
p6 <- show_sig_fit(optim_w, add.params = list(alpha = 0.1),
             samples = c("5115056-SRR8311885",
                         "5115615-SRR8311749",
                         "TCGA-KK-A7B4-01",
                         "TCGA-KK-A59X-01",
                         "TCGA-KC-A7FE-01"),
             plot_fun = "scatter",
             legend = "top") +
  scale_labels +
  ggplot2::guides(color=ggplot2::guide_legend(nrow=2,byrow=TRUE))
#pheatmap::pheatmap(optim_w)

library(patchwork)
p_all <- (p1 | p4) /
  p6 /
  (p2 | p3)
ggsave("figures_new/signature-fitting-and-bootstrap.pdf", plot = p_all, width = 11, height = 9)


# Copy number profile with method 'M' -------------------------------------

load("output/CNV.seqz.tally.M.RData")
load("output/Sig.CNV.seqz.M.RData")

p_m <- show_sig_profile(Sig.CNV.seqz.M,
                 mode = "copynumber", method = "M",
                 normalize = "column", style = "cosmic",
                 paint_axis_text = FALSE,
                 params = CNV.seqz.tally.M$parameters, y_expand = 1.5,
                 sig_names = paste("CN-Sig", 1:5))

p_w <- show_sig_profile(Sig.CNV.seqz.W$Signature[, 1, drop = FALSE],
                 mode = "copynumber", method = "W",
                 normalize = "feature", style = "cosmic",
                 sig_names = paste("CN-Sig", 1))

dir.create("figures_new")
ggsave(filename = "figures_new/signature_profile_m.pdf", plot = p_m, device = cairo_pdf,
       width = 10, height = 7)
ggsave(filename = "figures_new/signature_profile_w.pdf", plot = p_w, device = cairo_pdf,
       width = 12, height = 4)
