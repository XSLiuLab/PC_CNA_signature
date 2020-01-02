library(tidyverse)
library(sigminer)

df.facets = readRDS(file = "output/PRAD_Merge_Info_CNV_from_facets.RData")
df.seqz = readRDS(file = "output/PRAD_Merge_Info_CNV_from_sequenza.RData")

setdiff(colnames(df.facets), colnames(df.seqz))

cols_to_sigs <- c(paste0("CNV_Sig", 1:6), paste0("SNV_Sig", 1:6))
cols_to_sigs.seqz <- c(paste0("CNV_Sig", 1:5), paste0("SNV_Sig", 1:6))

# Exclude PSA
cols_to_features <- c(
  "IsMetastatic", "HasFusion",
  "Age", "Stage", "GleasonScore",
  "total_mutation", "n_driver", "Ti_fraction", "Tv_fraction",
  "n_of_cnv", "n_of_amp", "n_of_del", "n_of_vchr", "cna_burden",
  "MATH", "cluster", "purity", "ploidy"
)

cols_to_mutated_genes <- colnames(df.facets)[30:76]
cols_to_mutated_pathways <- colnames(df.facets)[77:87]

feature_type <- c(rep("ca", 2), rep("co", 16))

# FACETS ------------------------------------------------------------------

tidy_data.feature <- get_sig_feature_association(df.facets,
                                                 cols_to_sigs = cols_to_sigs,
                                                 cols_to_features = cols_to_features,
                                                 method_co = "pearson",
                                                 type = feature_type, verbose = TRUE) %>%
  get_tidy_association(p_adjust = TRUE)

show_sig_feature_corrplot(tidy_data.feature)

hist(df.facets$CNV_Sig1, breaks = 100)

tidy_data.gene <- get_sig_feature_association(df.facets,
                                              cols_to_sigs = cols_to_sigs,
                                              cols_to_features = cols_to_mutated_genes,
                                              type = "ca", verbose = TRUE) %>%
  get_tidy_association(p_adjust = TRUE)

show_sig_feature_corrplot(tidy_data.gene, ylab = "Mutated genes")

tidy_data.pathways <- get_sig_feature_association(df.facets,
                                                  cols_to_sigs = cols_to_sigs,
                                                  cols_to_features = cols_to_mutated_pathways,
                                                  type = "ca", verbose = TRUE) %>%
  get_tidy_association(p_adjust = TRUE)

show_sig_feature_corrplot(tidy_data.pathways, ylab = "Mutated pathways")

df.facets %>%
  dplyr::filter(!is.na(PSA)) %>%
  dplyr::pull(Study) %>%
  table()

df.psa = df.facets %>%
  dplyr::filter(!is.na(PSA)) %>%
  dplyr::select(contains("Sig", ignore.case = F), c("PSA", "Study", "CNV_ID"))

df.psa <- df.psa %>%
  dplyr:::group_by(Study) %>%
  tidyr::nest() %>%
  dplyr::summarise(
    data = purrr:::map(data, function(x) {
      get_sig_feature_association(x,
                                  cols_to_sigs = cols_to_sigs,
                                  cols_to_features = "PSA",
                                  method_co = "pearson",
                                  type = "co", verbose = TRUE) %>%
        get_tidy_association(p_adjust = TRUE)
    })
  ) %>%
  tidyr::unnest("data") %>%
  dplyr::filter(feature == "PSA") %>%
  dplyr::mutate(feature = Study)

show_sig_feature_corrplot(df.psa, p_val = 1, breaks_count = c(0, 50, 100, 150, 200), ylab = "PSA")

# Sequenza ----------------------------------------------------------------

tidy_data.seqz.feature <- get_sig_feature_association(df.seqz,
                                                      cols_to_sigs = cols_to_sigs.seqz,
                                                      cols_to_features = cols_to_features,
                                                      method_co = "pearson",
                                                      type = feature_type, verbose = TRUE) %>%
  get_tidy_association(p_adjust = TRUE)

show_sig_feature_corrplot(tidy_data.seqz.feature)

tidy_data.seqz.gene <- get_sig_feature_association(df.seqz,
                                                   cols_to_sigs = cols_to_sigs.seqz,
                                                   cols_to_features = cols_to_mutated_genes,
                                                   type = "ca", verbose = TRUE) %>%
  get_tidy_association(p_adjust = TRUE)

show_sig_feature_corrplot(tidy_data.seqz.gene, ylab = "Mutated genes")

tidy_data.seqz.pathways <- get_sig_feature_association(df.seqz,
                                                       cols_to_sigs = cols_to_sigs.seqz,
                                                       cols_to_features = cols_to_mutated_pathways,
                                                       type = "ca", verbose = TRUE) %>%
  get_tidy_association(p_adjust = TRUE)

show_sig_feature_corrplot(tidy_data.seqz.pathways, ylab = "Mutated pathways")


df.seqz.psa = df.seqz %>%
  dplyr::filter(!is.na(PSA)) %>%
  dplyr::select(contains("Sig", ignore.case = F), c("PSA", "Study", "CNV_ID"))

df.seqz.psa <- df.seqz.psa %>%
  dplyr:::group_by(Study) %>%
  tidyr::nest() %>%
  dplyr::summarise(
    data = purrr:::map(data, function(x) {
      get_sig_feature_association(x,
                                  cols_to_sigs = cols_to_sigs.seqz,
                                  cols_to_features = "PSA",
                                  method_co = "pearson",
                                  type = "co", verbose = TRUE) %>%
        get_tidy_association(p_adjust = TRUE)
    })
  ) %>%
  tidyr::unnest("data") %>%
  dplyr::filter(feature == "PSA") %>%
  dplyr::mutate(feature = Study)

show_sig_feature_corrplot(df.seqz.psa, p_val = 1, breaks_count = c(0, 50, 100, 150, 200), ylab = "PSA")
