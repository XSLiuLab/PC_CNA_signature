library(tidyverse)
library(sigminer)

load(file = "data/PRAD_Merge_Info.RData")
# Analyze association between signature and other features

cols_to_sigs <- c(paste0("CNV_Sig", 1:6), paste0("SNV_Sig", 1:3))
cols_to_features <- c(
  "IsMetastatic", "HasFusion", "Age", "Stage", "GleasonScore",
  "total_mutation", "n_driver", "Ti_fraction", "Tv_fraction", "MATH", "cluster",
  "n_of_cnv", "n_of_amp", "n_of_del", "cna_burden", "purity", "ploidy"
)
feature_type <- c(rep("ca", 2), rep("co", 15))

asso_data <- get_sig_feature_association(MergeInfo,
  cols_to_sigs = cols_to_sigs,
  cols_to_features = cols_to_features,
  method_co = "pearson",
  type = feature_type, verbose = TRUE
)

asso_tidy <- get_tidy_association(asso_data)

# Show corrplot
col3 <- colorRampPalette(c("red", "white", "blue"))
corrplot::corrplot(asso_data$corr_co$measure,
  col = rev(col3(20)),
  tl.col = "black",
  p.mat = asso_data$corr_co$p
)




show_sig_feature_corrplot(asso_tidy %>% dplyr::filter(p < 0.05), feature_list = c("IsMetastatic", "HasFusion"))
show_sig_feature_corrplot(asso_tidy %>% dplyr::filter(p < 0.05))
show_sig_feature_corrplot(asso_tidy)

# Show correlation network ------------------------------------------------

cols <- colnames(asso_data$corr_co$measure)
library(corrr)
res.cor <- correlate(
  MergeInfo[, cols] %>%
    dplyr::mutate_if(is.ordered, as.numeric)
)
# %>%
#   dplyr::mutate_if(is.numeric, dplyr::coalesce, 0)

# Talk this with package author
res.cor %>%
  network_plot(min_cor = 0.2, colours = rev(c("indianred2", "white", "skyblue1")))
# Remove variable causing plot failure
res.cor[c(-18, -21), c(-19, -22)] %>%
  network_plot(min_cor = 0.2, colours = rev(c("indianred2", "white", "skyblue1")))


# test_data = res.cor
# test_data$rowname = paste0("f", 1:24)
# colnames(test_data)[-1] = paste0("f", 1:24)
# #test_data[,-1] = test_data[, -1] * 0.999
# network_plot(test_data, min_cor = 0.2, colours = rev(c("indianred2", "white", "skyblue1")))
# tt %>%
#   network_plot(min_cor = 0.2, colours = rev(c("indianred2", "white", "skyblue1")))
# # Remove variable causing plot failure
# test_data[c(-18, -21), c(-19, -22)] %>%
#   network_plot(min_cor = 0.2, colours = rev(c("indianred2", "white", "skyblue1")))



# Select sample with most signature exposure and see their profile ---------
load(file = "output/CNV.prob.RData")
load("data/PRAD_CNV.RData")

Comp_df <- CNV.prob$nmf_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample") %>%
  as_tibble()

# Sig 4 enrich
show_cn_profile(
  data = CNV, chrs = paste0("chr", c(1:22, "X")), nrow = 3, ncol = 2, show_title = T,
  samples = MergeInfo %>%
    select(CNV_ID, cnv_group, cnv_enrich_sig, CNV_Sig1:CNV_Sig6) %>%
    filter(cnv_enrich_sig == "Sig4") %>%
    arrange(desc(CNV_Sig4)) %>%
    left_join(Comp_df %>% select(sample, starts_with("bpchrarm")), by = c("CNV_ID" = "sample")) %>%
    slice(1:6) %>% pull(CNV_ID)
)

# Sig4 enrich and bpchrarm2 high
show_cn_profile(
  data = CNV, chrs = paste0("chr", c(1:22, "X")), nrow = 3, ncol = 2, show_title = T,
  samples = MergeInfo %>%
    select(CNV_ID, cnv_group, cnv_enrich_sig, CNV_Sig1:CNV_Sig6) %>%
    filter(cnv_enrich_sig == "Sig4") %>%
    arrange(desc(CNV_Sig4)) %>%
    left_join(Comp_df %>% select(sample, starts_with("bpchrarm")), by = c("CNV_ID" = "sample")) %>%
    arrange(desc(bpchrarm2)) %>%
    slice(1:6) %>% pull(CNV_ID)
)


# Sig 1 enrich
show_cn_profile(
  data = CNV, chrs = paste0("chr", c(1:22, "X")), nrow = 3, ncol = 2, show_title = T,
  samples = MergeInfo %>%
    select(CNV_ID, cnv_group, cnv_enrich_sig, CNV_Sig1:CNV_Sig6) %>%
    filter(cnv_enrich_sig == "Sig1") %>%
    arrange(desc(CNV_Sig1)) %>%
    slice(1:6) %>% pull(CNV_ID)
)

# Sig 2 enrich
show_cn_profile(
  data = CNV, chrs = paste0("chr", c(1:22, "X")), nrow = 3, ncol = 2, show_title = T,
  samples = MergeInfo %>%
    select(CNV_ID, cnv_group, cnv_enrich_sig, CNV_Sig1:CNV_Sig6) %>%
    filter(cnv_enrich_sig == "Sig2") %>%
    arrange(desc(CNV_Sig2)) %>%
    slice(1:6) %>% pull(CNV_ID)
)

# Sig 3 enrich
show_cn_profile(
  data = CNV, chrs = paste0("chr", c(1:22, "X")), nrow = 3, ncol = 2, show_title = T,
  samples = MergeInfo %>%
    select(CNV_ID, cnv_group, cnv_enrich_sig, CNV_Sig1:CNV_Sig6) %>%
    filter(cnv_enrich_sig == "Sig3") %>%
    arrange(desc(CNV_Sig3)) %>%
    slice(1:6) %>% pull(CNV_ID)
)

# Sig 5 enrich
show_cn_profile(
  data = CNV, chrs = paste0("chr", c(1:22, "X")), nrow = 3, ncol = 2, show_title = T,
  samples = MergeInfo %>%
    select(CNV_ID, cnv_group, cnv_enrich_sig, CNV_Sig1:CNV_Sig6) %>%
    filter(cnv_enrich_sig == "Sig5") %>%
    arrange(desc(CNV_Sig5)) %>%
    slice(1:6) %>% pull(CNV_ID)
)

# Sig 6 enrich
show_cn_profile(
  data = CNV, chrs = paste0("chr", c(1:22, "X")), nrow = 3, ncol = 2, show_title = T,
  samples = MergeInfo %>%
    select(CNV_ID, cnv_group, cnv_enrich_sig, CNV_Sig1:CNV_Sig6) %>%
    filter(cnv_enrich_sig == "Sig6") %>%
    arrange(desc(CNV_Sig6)) %>%
    slice(1:6) %>% pull(CNV_ID)
)


# Sig 3 is very strange, extract the samples and check their bam files
special_samples <- MergeInfo %>%
  select(CNV_ID, cnv_group, cnv_enrich_sig, CNV_Sig1:CNV_Sig6) %>%
  filter(cnv_enrich_sig == "Sig3") %>%
  arrange(desc(CNV_Sig3)) %>%
  slice(1:6) %>%
  pull(CNV_ID)

Info <- readRDS("data/PRAD_CLINICAL.rds")
SpeInfo <- Info %>%
  filter(CNV_ID %in% special_samples) %>%
  select(CNV_ID, tumor_Run, normal_Run)
cat(SpeInfo$tumor_Run, SpeInfo$normal_Run)

SpeCNV <- CNV@data[sample %in% special_samples][segVal != 2]
readr::write_csv(SpeCNV, path = "data/CNV_for_sig3_top_enrich.csv")

o("/public/data/facet/dbGAP_PLUS_TCGA_PRAD_WES_CVAL150/915-1115081.pdf")

CVAL300 <- data.table::fread("data/CNV_from_dbGAP_PLUS_TCGA_WES_CVAL300.tsv")
CVAL300[sample %in% special_samples][modal_cn != 2]
# Show sig group comparison -----------------------------------------------

feature_type2 <- feature_type
feature_type2[4] <- "ca"
groups.cmp <- get_group_comparison(MergeInfo,
  col_group = "cnv_enrich_sig",
  cols_to_compare = cols_to_features,
  type = feature_type2, verbose = TRUE
)
p.cmps <- show_group_comparison(
  group_comparison = groups.cmp,
  xlab = NA,
  method = "anova",
  legend_position_ca = "right",
  hjust = 1,
  text_angle_x = 60
)
p.cmps$co_comb
p.cmps$ca_comb

# For mutational signatures
groups.cmp.mt <- get_group_comparison(MergeInfo,
  col_group = "snv_enrich_sig",
  cols_to_compare = cols_to_features,
  type = feature_type2, verbose = TRUE
)
p.cmps.mt <- show_group_comparison(
  group_comparison = groups.cmp.mt,
  xlab = NA,
  method = "anova",
  legend_position_ca = "right",
  hjust = 1,
  text_angle_x = 60
)
p.cmps.mt$co_comb
p.cmps.mt$ca_comb


# Show signature group maps -----------------------------------------------

show_group_mapping(MergeInfo,
  sig_col = "cnv_enrich_sig",
  map_cols = c("GleasonScore"),
  include_sig = TRUE
)

show_group_mapping(MergeInfo[cnv_enrich_sig %in% paste0("Sig", 1:6)],
  sig_col = "cnv_enrich_sig",
  map_cols = c("snv_enrich_sig", "sample_type", "GleasonScore", "Stage"),
  include_sig = TRUE
)

show_group_mapping(MergeInfo[cnv_enrich_sig %in% paste0("Sig", 1:6)],
  sig_col = "cnv_enrich_sig",
  map_cols = c("sample_type", "snv_enrich_sig"),
  include_sig = TRUE
)
