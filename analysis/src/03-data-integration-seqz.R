# Integrate all informaiton to sample level
library(tidyverse)
library(sigminer)
library(maftools)

# Loading clinical related data -------------------------------------------

Info <- readRDS("data/PRAD_CLINICAL.rds")
# Purity and ploidy info from facets
PurityInfo <- read_tsv("data/PRAD_Purity_and_Ploidy_Sequenza.tsv")

# Processing CNV data -----------------------------------------------------

load("output/CNV.seqz.RData")
#load("output/Sig.CNV.seqz.W.RData")
load("output/Sig.CNV.seqz.W.5.RData")
CNV <- CNV.seqz
Sig.CNV <- Sig.CNV.seqz.W.5
rm(Sig.CNV.seqz.W.5, CNV.seqz)

# CNV
CNVGroupInfo <- get_groups(Sig.CNV, method = "consensus", match_consensus = TRUE)
CNVInfo <- CNV@summary.per.sample
CNVExposureInfo <- get_sig_exposure(Sig.CNV)

# Processing mutation data ------------------------------------------------

load(file = "output/PRAD_TCGA_plus_dbGap_Maf.RData")
load(file = "output/Sig.PRAD_TCGA_plus_dbGap_rm_hyper.RData")

TMBInfo <- getSampleSummary(Maf)[, .(Tumor_Sample_Barcode, total)]

load(file = "output/PRAD_driver_info.RData")
load(file = "output/PRAD_heter_info.RData")

SNVGroupInfo <- get_groups(Sig.SNV, method = "consensus", match_consensus = TRUE)
SNVExposureInfo <- get_sig_exposure(Sig.SNV)

# Processing gene and pathway mutation ------------------------------------

load(file = "output/PRAD_gene_and_pathway_mutation.RData")

# Merge data --------------------------------------------------------------
Info <- Info %>%
  dplyr::mutate(
    CNV_ID = dplyr::case_when(
      !startsWith(tumor_Run, "TCGA") & !is.na(tumor_Run) ~ paste(subject_id, tumor_Run, sep = "-"),
      startsWith(tumor_Run, "TCGA") & !is.na(tumor_Run) ~ tumor_Run,
      TRUE ~ NA_character_
    )
  )

PurityInfo
colnames(CNVGroupInfo) <- c("sample", "cnv_group", "cnv_weight", "cnv_enrich_sig")
CNVInfo
colnames(CNVExposureInfo) <- c("sample", paste0("CNV_", colnames(CNVExposureInfo)[-1]))
colnames(TMBInfo) <- c("sample", "total_mutation")
TMBInfo[, sample := as.character(sample)][, sample := ifelse(startsWith(sample, "TCGA"),
  substr(sample, 1, 15),
  sample
)]
DriverDF[, sample := as.character(sample)][, sample := ifelse(startsWith(sample, "TCGA"),
  substr(sample, 1, 15),
  sample
)]
colnames(SNVGroupInfo) <- c("sample", "snv_group", "snv_weight", "snv_enrich_sig")
SNVGroupInfo[, sample := as.character(sample)][, sample := ifelse(startsWith(sample, "TCGA"),
  substr(sample, 1, 15),
  sample
)]
colnames(SNVExposureInfo) <- c("sample", paste0("SNV_", colnames(SNVExposureInfo)[-1]))
SNVExposureInfo[, sample := as.character(sample)][, sample := ifelse(startsWith(sample, "TCGA"),
  substr(sample, 1, 15),
  sample
)]
colnames(TitvInfo) <- c("sample", "Ti_fraction", "Tv_fraction")
TitvInfo[, sample := as.character(sample)][, sample := ifelse(startsWith(sample, "TCGA"),
  substr(sample, 1, 15),
  sample
)]
colnames(MathDF) <- c("sample", "MATH")
MathDF[, sample := as.character(sample)][, sample := ifelse(startsWith(sample, "TCGA"),
  substr(sample, 1, 15),
  sample
)]
colnames(ClusterDF) <- c("sample", "cluster")
ClusterDF[, sample := as.character(sample)][, sample := ifelse(startsWith(sample, "TCGA"),
  substr(sample, 1, 15),
  sample
)]

data.table::setDT(summary_mutation)
data.table::setDT(summary_pathway)

colnames(summary_mutation)[1] <- "sample"
colnames(summary_pathway)[1] <- "sample"
summary_mutation[, sample := as.character(sample)][, sample := ifelse(startsWith(sample, "TCGA"),
  substr(sample, 1, 15),
  sample
)]
summary_pathway[, sample := as.character(sample)][, sample := ifelse(startsWith(sample, "TCGA"),
  substr(sample, 1, 15),
  sample
)]


MergeInfo <- Info %>%
  left_join(CNVInfo, by = c("CNV_ID" = "sample")) %>%
  left_join(CNVGroupInfo, by = c("CNV_ID" = "sample")) %>%
  left_join(CNVExposureInfo, by = c("CNV_ID" = "sample")) %>%
  left_join(PurityInfo, by = c("CNV_ID" = "sample")) %>%
  left_join(summary_mutation, by = c("tumor_Run" = "sample")) %>%
  left_join(summary_pathway, by = c("tumor_Run" = "sample")) %>%
  left_join(TMBInfo, by = c("tumor_Run" = "sample")) %>%
  left_join(DriverDF, by = c("tumor_Run" = "sample")) %>%
  dplyr::mutate(
    n_driver = ifelse(!is.na(n_driver), n_driver, 0)
  ) %>%
  left_join(TitvInfo, by = c("tumor_Run" = "sample")) %>%
  left_join(MathDF, by = c("tumor_Run" = "sample")) %>%
  left_join(ClusterDF, by = c("tumor_Run" = "sample")) %>%
  left_join(SNVGroupInfo, by = c("tumor_Run" = "sample")) %>%
  left_join(SNVExposureInfo, by = c("tumor_Run" = "sample")) %>%
  mutate(
    Stage = factor(Stage, ordered = TRUE),
    Fusion = ifelse(Fusion == "Negative", "No", "Yes"),
    sample_type = ifelse(sample_type == "Unknown", NA_character_, sample_type),
    HasFusion = Fusion,
    HasFusion = ifelse(HasFusion == "Yes", TRUE, FALSE),
    IsMetastatic = ifelse(sample_type == "Metastatic", TRUE, FALSE)
  )

summary(MergeInfo)
saveRDS(MergeInfo, file = "output/PRAD_Merge_Info_CNV_from_sequenza.RData")
