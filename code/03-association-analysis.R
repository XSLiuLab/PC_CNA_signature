# Analyze association between signature and other features

cols_to_sigs = c(paste0('CNV_Sig',1:6), paste0('SNV_Sig',1:3))
cols_to_features = c('IsMetastatic', 'HasFusion', 'Age', 'Stage', 'GleasonScore',
                     'total_mutation', "n_driver", "Ti_fraction", "Tv_fraction", "MATH", "cluster",
                     'n_of_cnv', 'n_of_amp', 'n_of_del', 'cna_burden', 'purity', 'ploidy')
feature_type = c(rep('ca', 2), rep('co', 15))

asso_data = get_sig_feature_association(MergeInfo,
                                        cols_to_sigs = cols_to_sigs,
                                        cols_to_features = cols_to_features,
                                        method_co = "pearson",
                                        type = feature_type, verbose = TRUE)

asso_tidy = get_tidy_association(asso_data)

# Show corrplot
col3 <- colorRampPalette(c("red", "white", "blue"))
corrplot::corrplot(asso_data$corr_co$measure,
                   col = rev(col3(20)),
                   tl.col = "black",
                   p.mat = asso_data$corr_co$p
)


show_sig_feature_corrplot(asso_tidy, features_list = c("IsMetastatic", "HasFusion"))



# Show correlation network ------------------------------------------------

cols = colnames(asso_data$corr_co$measure)
library(corrr)
res.cor = correlate(
  MergeInfo[, cols] %>%
    dplyr::mutate_if(is.ordered, as.numeric)
) %>%
  dplyr::mutate_if(is.numeric, dplyr::coalesce, 0)

# Remove variable causing plot failure
res.cor[-19, -20] %>%
  network_plot(min_cor = 0.2, colours = rev(c("indianred2", "white", "skyblue1")))

# Select sample with most signature exposure and see their profile

Comp_df = CNV.prob$nmf_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample") %>%
  as_tibble()

# Sig 4 enrich
show_cn_profile(data = CNV, chrs = paste0("chr", c(1:22, "X")), nrow = 3, ncol = 2, show_title = T,
                samples = MergeInfo %>%
                  select(CNV_ID, cnv_group, cnv_enrich_sig, CNV_Sig1:CNV_Sig6) %>%
                  filter(cnv_enrich_sig == 'Sig4') %>%
                  arrange(desc(CNV_Sig4)) %>%
                  left_join(Comp_df %>% select(sample, starts_with("bpchrarm")), by = c("CNV_ID"="sample")) %>%
                  slice(1:6) %>% pull(CNV_ID))

# Sig4 enrich and bpchrarm2 high
show_cn_profile(data = CNV, chrs = paste0("chr", c(1:22, "X")), nrow = 3, ncol = 2, show_title = T,
                samples =MergeInfo %>%
                  select(CNV_ID, cnv_group, cnv_enrich_sig, CNV_Sig1:CNV_Sig6) %>%
                  filter(cnv_enrich_sig == 'Sig4') %>%
                  arrange(desc(CNV_Sig4)) %>%
                  left_join(Comp_df %>% select(sample, starts_with("bpchrarm")), by = c("CNV_ID"="sample")) %>%
                  arrange(desc(bpchrarm2)) %>%
                  slice(1:6) %>% pull(CNV_ID))


# Sig 1 enrich
show_cn_profile(data = CNV, chrs = paste0("chr", c(1:22, "X")), nrow = 3, ncol = 2, show_title = T,
                samples = MergeInfo %>%
                  select(CNV_ID, cnv_group, cnv_enrich_sig, CNV_Sig1:CNV_Sig6) %>%
                  filter(cnv_enrich_sig == 'Sig1') %>%
                  arrange(desc(CNV_Sig1)) %>%
                  slice(1:6) %>% pull(CNV_ID))

# Sig 2 enrich
show_cn_profile(data = CNV, chrs = paste0("chr", c(1:22, "X")), nrow = 3, ncol = 2, show_title = T,
                samples = MergeInfo %>%
                  select(CNV_ID, cnv_group, cnv_enrich_sig, CNV_Sig1:CNV_Sig6) %>%
                  filter(cnv_enrich_sig == 'Sig2') %>%
                  arrange(desc(CNV_Sig2)) %>%
                  slice(1:6) %>% pull(CNV_ID))

# Sig 3 enrich
show_cn_profile(data = CNV, chrs = paste0("chr", c(1:22, "X")), nrow = 3, ncol = 2, show_title = T,
                samples = MergeInfo %>%
                  select(CNV_ID, cnv_group, cnv_enrich_sig, CNV_Sig1:CNV_Sig6) %>%
                  filter(cnv_enrich_sig == 'Sig3') %>%
                  arrange(desc(CNV_Sig3)) %>%
                  slice(1:6) %>% pull(CNV_ID))

# Sig 5 enrich
show_cn_profile(data = CNV, chrs = paste0("chr", c(1:22, "X")), nrow = 3, ncol = 2, show_title = T,
                samples = MergeInfo %>%
                  select(CNV_ID, cnv_group, cnv_enrich_sig, CNV_Sig1:CNV_Sig6) %>%
                  filter(cnv_enrich_sig == 'Sig5') %>%
                  arrange(desc(CNV_Sig5)) %>%
                  slice(1:6) %>% pull(CNV_ID))

# Sig 6 enrich
show_cn_profile(data = CNV, chrs = paste0("chr", c(1:22, "X")), nrow = 3, ncol = 2, show_title = T,
                samples = MergeInfo %>%
                  select(CNV_ID, cnv_group, cnv_enrich_sig, CNV_Sig1:CNV_Sig6) %>%
                  filter(cnv_enrich_sig == 'Sig6') %>%
                  arrange(desc(CNV_Sig6)) %>%
                  slice(1:6) %>% pull(CNV_ID))
