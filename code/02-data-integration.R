# Integrate all informaiton to sample level
library(tidyverse)
library(sigminer)

load("data/PRAD_CNV.RData")
Info = readRDS("data/PRAD_CLINICAL.rds")
PurityInfo = read_tsv("data/PRAD_Purity_and_Ploidy_CVAL150.tsv")


CNV.Sig = readRDS("output/NMF_copynumber_signature.prob.rds")
#CNV.Sig.count = readRDS("output/NMF_copynumber_signature.count.rds")

# # Observe profile
load(file = "output/CNV.prob.RData")
load(file = "output/CNV.count.RData")
#
# show_sig_profile(CNV.Sig, params = CNV.prob$parameters, y_expand = 1.5)
# show_sig_profile(CNV.Sig, params = CNV.prob$parameters, y_expand = 1.5, normalize = 'column')
# show_sig_profile(CNV.Sig.count, params = CNV.count$parameters, y_expand = 1.5)
# show_sig_profile(CNV.Sig.count, params = CNV.count$parameters, y_expand = 1.5, normalize = 'column')


GroupInfo = get_groups(CNV.Sig)
CNVInfo = CNV@summary.per.sample
ExposureInfo = get_sig_exposure(CNV.Sig)

MergeInfo = Info %>%
  left_join(CNVInfo, by = c('CNV_ID'='sample')) %>%
  left_join(GroupInfo, by = c('CNV_ID' = 'sample')) %>%
  left_join(ExposureInfo, by = c('CNV_ID' = 'sample')) %>%
  select(Study, sample_type, PatientID:Sig6) %>%
  mutate(Stage = factor(Stage, ordered = TRUE),
         Fusion = ifelse(Fusion=='Negative', 'No', 'Yes'),
         sample_type = ifelse(sample_type == "Unknown", NA_character_, sample_type),
         HasFusion = Fusion,
         HasFusion = ifelse(HasFusion == 'Yes', TRUE, FALSE),
         IsMetastatic = ifelse(sample_type == 'Metastatic', TRUE, FALSE)) %>%
  left_join(PurityInfo, by = c('CNV_ID' = 'sample'))

summary(MergeInfo)

cols_to_sigs = paste0('Sig',1:6)
cols_to_features = c('IsMetastatic', 'HasFusion', 'Age', 'Stage', 'GleasonScore', 'n_of_cnv', 'n_of_amp', 'n_of_del', 'cna_burden', 'purity', 'ploidy')
feature_type = c(rep('ca', 2), rep('co', 9))

asso_data = get_sig_feature_association(MergeInfo,
                                        cols_to_sigs = cols_to_sigs,
                                        cols_to_features = cols_to_features,
                                        method_co = "spearman",
                                        type = feature_type, verbose = TRUE)
# Take a check
cor(MergeInfo$Sig2, MergeInfo$Sig1, use = 'pairwise.complete.obs', method = 'spearman')

asso_tidy = get_tidy_association(asso_data)

# Show corrplot
col3 <- colorRampPalette(c("red", "white", "blue"))
corrplot::corrplot(asso_data$corr_co$measure,
                   col = rev(col3(20)),
                   tl.col = "black",
                   p.mat = asso_data$corr_co$p
)


# Select sample with most signature exposure and see their profile

Comp_df = CNV.prob$nmf_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample") %>%
  as_tibble()

MergeInfo %>%
  select(CNV_ID, group, enrich_sig, Sig1:Sig6) %>%
  filter(enrich_sig == 'Sig5') %>%
  arrange(desc(Sig5)) %>%
  left_join(Comp_df %>% select(sample, starts_with("bpchrarm")), by = c("CNV_ID"="sample"))

o("/public/data/facet/dbGAP_PLUS_TCGA_PRAD_WES_CVAL150/TCGA-EJ-A46B.pdf")

MergeInfo %>%
  select(CNV_ID, group, enrich_sig, Sig1:Sig6) %>%
  filter(enrich_sig == 'Sig5') %>%
  arrange(desc(Sig5)) %>%
  left_join(Comp_df %>% select(sample, starts_with("bpchrarm")), by = c("CNV_ID"="sample")) %>%
  arrange(desc(bpchrarm2))

o("/public/data/facet/dbGAP_PLUS_TCGA_PRAD_WES_CVAL150/915-6115124.pdf")
o("/public/data/facet/dbGAP_PLUS_TCGA_PRAD_WES_CVAL150/909-21020.pdf")
o("/public/data/facet/dbGAP_PLUS_TCGA_PRAD_WES_CVAL150/554-4.pdf")
o("/public/data/facet/dbGAP_PLUS_TCGA_PRAD_WES_CVAL150/TCGA-EJ-5515.pdf")
o("/public/data/facet/dbGAP_PLUS_TCGA_PRAD_WES_CVAL150/915-5115076.pdf")

