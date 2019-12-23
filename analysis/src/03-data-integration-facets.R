# Integrate all informaiton to sample level
library(tidyverse)
library(sigminer)
library(maftools)

# Loading clinical related data -------------------------------------------

Info <- readRDS("data/PRAD_CLINICAL.rds")
# Purity and ploidy info from sequenza
PurityInfo <- read_tsv("data/PRAD_Purity_and_Ploidy_CVAL150.tsv")

# Processing CNV data -----------------------------------------------------

load("output/CNV.facets.RData")
load("output/Sig.CNV.facets.W.RData")
CNV <- CNV.facets
Sig.CNV <- Sig.CNV.facets.W
rm(Sig.CNV.facets.W, CNV.facets)

# CNV
CNVGroupInfo <- get_groups(Sig.CNV, method = "consensus", match_consensus = TRUE)
CNVInfo <- CNV@summary.per.sample
CNVExposureInfo <- get_sig_exposure(Sig.CNV)

# Processing mutation data ------------------------------------------------

load(file = "output/PRAD_TCGA_plus_dbGap_Maf.RData")
load(file = "output/Sig.PRAD_TCGA_plus_dbGap_rm_hyper.RData")

get_sig_similarity(Sig.SNV)

TMBInfo <- getSampleSummary(Maf)[, .(Tumor_Sample_Barcode, total)]

# Using oncodrive method
#
# DriverInfo = oncodrive(Maf)
# driver_genes = DriverInfo[fdr < 0.05 & pval < 0.01]$Hugo_Symbol

# Using MutSig method
DriverInfo <- data.table::fread("data/PRAD.sig_genes.txt")
driver_genes <- DriverInfo[p < 0.01 & q < 0.05]$gene

DriverDF <- map_df(genesToBarcodes(Maf, genes = driver_genes), function(x) {
  dplyr::tibble(sample = x$Tumor_Sample_Barcode)
}) %>%
  count(sample) %>%
  rename(n_driver = n)
data.table::setDT(DriverDF)

save(DriverInfo, driver_genes, DriverDF, file = "output/PRAD_driver_info.RData")

TitvInfo <- titv(maf = Maf, plot = FALSE, useSyn = TRUE)$TiTv.fractions
MathInfo <- inferHeterogeneity(Maf, TitvInfo$Tumor_Sample_Barcode, useSyn = TRUE)
MathDF <- MathInfo$clusterData[, list(MATH = mean(MATH, na.rm = TRUE)), by = Tumor_Sample_Barcode]
ClusterDF <- MathInfo$clusterMeans[, list(cluster = as.integer(cluster), Tumor_Sample_Barcode)][, list(cluster = max(cluster, na.rm = TRUE)), by = Tumor_Sample_Barcode]

save(TitvInfo, MathInfo, MathDF, ClusterDF, file = "output/PRAD_heter_info.RData")

SNVGroupInfo <- get_groups(Sig.SNV, method = "consensus", match_consensus = TRUE)
SNVExposureInfo <- get_sig_exposure(Sig.SNV)

# Processing gene and pathway mutation ------------------------------------

get_mutation_status <- function(x, genes) {
  as_tibble(
    setNames(
      lapply(genes, function(gene) any(grepl(paste0("^", gene, "$"), x))),
      genes
    )
  )
}

summary_mutation <- Maf@data %>%
  group_by(Tumor_Sample_Barcode) %>%
  do(get_mutation_status(.$Hugo_Symbol, driver_genes))


# Only functional variants included here
# get mutation status of pathways
get_pathway_mut_status <- function(x, paths) {
  as_tibble(
    setNames(
      lapply(paths, function(path) any(grepl(paste0("^", path, "$"), x))),
      paths
    )
  )
}

pathdb <- system.file("extdata", "oncogenic_sig_patwhays.tsv", package = "maftools")
pathdb <- data.table::fread(input = pathdb)

HR <- data.table::data.table(
  Pathway = "HR_pathway",
  Gene = c("BARD1", "PALB2", "BRCA1", "BRCA2", "ATR", "BLM", "ATM", "NBN", "MRE11")
)

pathdb <- rbind(pathdb, HR, fill = TRUE)

path_list <- unique(pathdb$Pathway)
mut_data <- dplyr::left_join(Maf@data, pathdb, by = c("Hugo_Symbol" = "Gene"))

summary_pathway <- mut_data %>%
  dplyr::group_by(Tumor_Sample_Barcode) %>%
  dplyr::do(get_pathway_mut_status(.$Pathway, path_list))
# Make TP53 pathway differ from TP53 gene
colnames(summary_pathway)[colnames(summary_pathway) == "TP53"] <- "TP53_pathway"

save(summary_mutation, summary_pathway, file = "output/PRAD_gene_and_pathway_mutation.RData")

# Merge data --------------------------------------------------------------
Info
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
saveRDS(MergeInfo, file = "output/PRAD_Merge_Info_CNV_from_facets.RData")
