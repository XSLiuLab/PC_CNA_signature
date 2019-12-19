# Integrate all informaiton to sample level
library(tidyverse)
library(sigminer)
library(maftools)

# Loading clinical related data -------------------------------------------

Info = readRDS("data/PRAD_CLINICAL.rds")
PurityInfo = read_tsv("data/PRAD_Purity_and_Ploidy_CVAL150.tsv")

# Processing CNV data -----------------------------------------------------

load("data/PRAD_CNV.RData")
CNV.Sig = readRDS("output/NMF_copynumber_signature.prob.rds")
#CNV.Sig.count = readRDS("output/NMF_copynumber_signature.count.rds")

# Observe profile
load(file = "output/CNV.prob.RData")
#load(file = "output/CNV.count.RData")
#
show_cn_distribution(CNV)
show_cn_features(CNV.prob$features)
show_cn_components(CNV.prob$parameters, show_weights = F)

show_sig_profile(CNV.Sig, params = CNV.prob$parameters, y_expand = 3, show_cv = TRUE, params_label_angle = 80)
show_sig_profile(CNV.Sig, params = CNV.prob$parameters, y_expand = 3, show_cv = TRUE, params_label_angle = 80,
                 normalize = 'column')
# show_sig_profile(CNV.Sig.count, params = CNV.count$parameters, y_expand = 1.5)
# show_sig_profile(CNV.Sig.count, params = CNV.count$parameters, y_expand = 1.5, normalize = 'column')

# CNV
CNVGroupInfo = get_groups(CNV.Sig)
CNVInfo = CNV@summary.per.sample
CNVExposureInfo = get_sig_exposure(CNV.Sig)

# Processing mutation data ------------------------------------------------

load("data/PRAD_Maf.RData")
SNV.Sig = readRDS(file = "output/NMF_snv_signature.rds")

show_sig_profile(SNV.Sig, mode = "mutation", normalize = "row")
get_sig_similarity(SNV.Sig)

TMBInfo = getSampleSummary(Maf)[, .(Tumor_Sample_Barcode, total)]

DriverInfo = oncodrive(Maf)
driver_genes = DriverInfo[fdr < 0.05 & pval < 0.01]$Hugo_Symbol
DriverDF = map_df(genesToBarcodes(Maf, genes = driver_genes), function(x) {
  dplyr::tibble(sample=x$Tumor_Sample_Barcode)
}) %>%
  count(sample) %>%
  rename(n_driver = n)
data.table::setDT(DriverDF)

TitvInfo = titv(maf = Maf, plot = FALSE, useSyn = TRUE)$TiTv.fractions
MathInfo = inferHeterogeneity(Maf, TitvInfo$Tumor_Sample_Barcode, useSyn = TRUE)
MathDF = MathInfo$clusterData[, list(MATH=mean(MATH, na.rm=TRUE)), by=Tumor_Sample_Barcode]
ClusterDF = MathInfo$clusterMeans[, list(cluster = as.integer(cluster),Tumor_Sample_Barcode)][, list(cluster=max(cluster, na.rm = TRUE)), by=Tumor_Sample_Barcode]

SNVGroupInfo = get_groups(SNV.Sig)
SNVExposureInfo = get_sig_exposure(SNV.Sig)



# Processing gene and pathway mutation ------------------------------------



# Merge data --------------------------------------------------------------
Info
PurityInfo
colnames(CNVGroupInfo) = c("sample", "cnv_group", "cnv_weight", "cnv_enrich_sig")
CNVInfo
colnames(CNVExposureInfo) = c("sample", paste0("CNV_Sig", 1:6))
colnames(TMBInfo) = c("sample", "total_mutation")
TMBInfo[, sample:=as.character(sample)][, sample:=ifelse(startsWith(sample, "TCGA"),
                         substr(sample, 1, 15),
                         sample)]
DriverDF[, sample:=as.character(sample)][, sample:=ifelse(startsWith(sample, "TCGA"),
                                                          substr(sample, 1, 15),
                                                          sample)]
colnames(SNVGroupInfo) = c("sample", "snv_group", "snv_weight", "snv_enrich_sig")
SNVGroupInfo[, sample:=as.character(sample)][, sample:=ifelse(startsWith(sample, "TCGA"),
                                                              substr(sample, 1, 15),
                                                              sample)]
colnames(SNVExposureInfo) = c("sample", paste0("SNV_Sig", 1:3))
SNVExposureInfo[, sample:=as.character(sample)][, sample:=ifelse(startsWith(sample, "TCGA"),
                                                                 substr(sample, 1, 15),
                                                                 sample)]
colnames(TitvInfo) = c("sample", "Ti_fraction", "Tv_fraction")
TitvInfo[, sample:=as.character(sample)][, sample:=ifelse(startsWith(sample, "TCGA"),
                                                          substr(sample, 1, 15),
                                                          sample)]
colnames(MathDF) = c("sample", "MATH")
MathDF[, sample:=as.character(sample)][, sample:=ifelse(startsWith(sample, "TCGA"),
                                                  substr(sample, 1, 15),
                                                  sample)]
colnames(ClusterDF) = c("sample", "cluster")
ClusterDF[, sample:=as.character(sample)][, sample:=ifelse(startsWith(sample, "TCGA"),
                                                           substr(sample, 1, 15),
                                                           sample)]

MergeInfo = Info %>%
  left_join(CNVInfo, by = c('CNV_ID'='sample')) %>%
  left_join(CNVGroupInfo, by = c('CNV_ID' = 'sample')) %>%
  left_join(CNVExposureInfo, by = c('CNV_ID' = 'sample')) %>%
  left_join(TMBInfo, by = c('tumor_Run'='sample')) %>%
  left_join(DriverDF, by = c('tumor_Run'='sample')) %>%
  dplyr::mutate(
    n_driver = ifelse(!is.na(n_driver), n_driver, 0)
  ) %>%
  left_join(TitvInfo, by = c('tumor_Run'='sample')) %>%
  left_join(MathDF, by = c('tumor_Run'='sample')) %>%
  left_join(ClusterDF, by = c('tumor_Run'='sample')) %>%
  left_join(SNVGroupInfo, by = c('tumor_Run'='sample')) %>%
  left_join(SNVExposureInfo, by = c('tumor_Run'='sample')) %>%
  select(Study, sample_type, PatientID:SNV_Sig3) %>%
  mutate(Stage = factor(Stage, ordered = TRUE),
         Fusion = ifelse(Fusion=='Negative', 'No', 'Yes'),
         sample_type = ifelse(sample_type == "Unknown", NA_character_, sample_type),
         HasFusion = Fusion,
         HasFusion = ifelse(HasFusion == 'Yes', TRUE, FALSE),
         IsMetastatic = ifelse(sample_type == 'Metastatic', TRUE, FALSE)) %>%
  left_join(PurityInfo, by = c('CNV_ID' = 'sample'))

summary(MergeInfo)

save(MergeInfo, file = "data/PRAD_Merge_Info.RData")



# Show signature exposure -------------------------------------------------

show_sig_exposure(CNV.Sig, rm_space = TRUE)
show_sig_exposure(SNV.Sig, rm_space = TRUE)
