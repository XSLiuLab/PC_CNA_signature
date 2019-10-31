# Integrate all informaiton to sample level
library(tidyverse)
library(sigminer)

load("data/PRAD_CNV.RData")
Info = readRDS("data/PRAD_CLINICAL.rds")
CNV.Sig = readRDS("output/NMF_copynumber_signature.prob.rds")
CNV.Sig.count = readRDS("output/NMF_copynumber_signature.count.rds")

# Observe profile
load(file = "output/CNV.prob.RData")
load(file = "output/CNV.count.RData")

show_sig_profile(CNV.Sig, params = CNV.prob$parameters, y_expand = 1.5)
show_sig_profile(CNV.Sig, params = CNV.prob$parameters, y_expand = 1.5, normalize = 'column')
show_sig_profile(CNV.Sig.count, params = CNV.prob$parameters, y_expand = 1.5)
show_sig_profile(CNV.Sig.count, params = CNV.prob$parameters, y_expand = 1.5, normalize = 'column')


GroupInfo = get_groups(CNV.Sig, match_consensus = TRUE)
CNVInfo = CNV@summary.per.sample
ExposureInfo = get_sig_exposure(CNV.Sig)

MergeInfo = Info %>%
  left_join(CNVInfo, by = c('CNV_ID'='sample')) %>%
  left_join(GroupInfo, by = c('CNV_ID' = 'sample')) %>%
  left_join(ExposureInfo, by = c('CNV_ID' = 'sample')) %>%
  select(Study, sample_type, PatientID:Sig5) %>%
  mutate(Stage = factor(Stage, ordered = TRUE),
         Fusion = ifelse(Fusion=='Negative', 'No', 'Yes'),
         sample_type = ifelse(sample_type == "Unknown", NA_character_, sample_type),
         HasFusion = Fusion,
         HasFusion = ifelse(HasFusion == 'Yes', TRUE, FALSE),
         IsMetastatic = ifelse(sample_type == 'Metastatic', TRUE, FALSE))

summary(MergeInfo)

cols_to_sigs = paste0('Sig',1:5)
cols_to_features = c('IsMetastatic', 'HasFusion', 'Age', 'Stage', 'GleasonScore', 'n_of_cnv', 'n_of_amp', 'n_of_del', 'cna_burden')
feature_type = c(rep('ca', 2), rep('co', 7))

asso_data = get_sig_feature_association(MergeInfo,
                                        cols_to_sigs = cols_to_sigs,
                                        cols_to_features = cols_to_features,
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
