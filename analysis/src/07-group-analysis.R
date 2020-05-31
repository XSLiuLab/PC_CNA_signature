## Old code

library(sigminer)

load("output/CNV.seqz.RData")
load("output/Sig.CNV.seqz.W.5.RData")
load("output/CNV.seqz.derive.W.RData")
# Copy number signature dominant samples ----------------------------------

cnv_group <- get_groups(Sig.CNV.seqz.W.5)
cnv_expo <- get_sig_exposure(Sig.CNV.seqz.W.5, type = "relative")

df <- dplyr::left_join(cnv_group, cnv_expo)

show_cn_profile(
  data = CNV.seqz, nrow = 3, ncol = 2, show_title = T,
  samples = df %>%
    dplyr::filter(enrich_sig == "Sig1") %>%
    dplyr::arrange(dplyr::desc(Sig1)) %>%
    dplyr::slice(1:6) %>% dplyr::pull(sample)
)

show_cn_profile(
  data = CNV.seqz, nrow = 3, ncol = 2, show_title = T,
  samples = df %>%
    dplyr::filter(enrich_sig == "Sig2") %>%
    dplyr::arrange(dplyr::desc(Sig2)) %>%
    dplyr::slice(1:6) %>% dplyr::pull(sample)
)

show_cn_profile(
  data = CNV.seqz, nrow = 3, ncol = 2, show_title = T,
  samples = df %>%
    dplyr::filter(enrich_sig == "Sig3") %>%
    dplyr::arrange(dplyr::desc(Sig3)) %>%
    dplyr::slice(1:6) %>% dplyr::pull(sample)
)

show_cn_profile(
  data = CNV.seqz, nrow = 3, ncol = 2, show_title = T,
  samples = df %>%
    dplyr::filter(enrich_sig == "Sig4") %>%
    dplyr::arrange(dplyr::desc(Sig4)) %>%
    dplyr::slice(1:6) %>% dplyr::pull(sample)
)

show_cn_profile(
  data = CNV.seqz, nrow = 3, ncol = 2, show_title = T,
  samples = df %>%
    dplyr::filter(enrich_sig == "Sig5") %>%
    dplyr::arrange(dplyr::desc(Sig5)) %>%
    dplyr::slice(1:6) %>% dplyr::pull(sample)
)



show_cn_profile(
  data = CNV.seqz, nrow = 3, ncol = 2, show_title = T,
  samples = df %>%
    dplyr::filter(enrich_sig == "Sig5") %>%
    dplyr::left_join(as.data.frame(CNV.seqz.derive.W$nmf_matrix) %>%
                tibble::rownames_to_column("sample") %>%
                dplyr::select(sample, `BPArm[2]`, `BPArm[1]`)) %>%
    dplyr::mutate(BK = `BPArm[2]` + `BPArm[1]`) %>%
    dplyr::arrange(desc(BK)) %>%
    dplyr::slice(1:6) %>% dplyr::pull(sample)
)
