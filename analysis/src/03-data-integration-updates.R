# Integrate all informaiton to sample level
library(tidyverse)
library(sigminer)
library(maftools)

# Loading clinical related data -------------------------------------------

Info <- readRDS("data/PRAD_CLINICAL.rds")
# Purity and ploidy info from sequenza
PurityInfo <- read_tsv("data/PRAD_Purity_and_Ploidy_Sequenza.tsv")

# Processing CNV data -----------------------------------------------------

load("output/CNV.seqz.RData")
load("output/Sig.CNV.seqz.W.RData")
CNV <- CNV.seqz
Sig.CNV <- Sig.CNV.seqz.W
rm(Sig.CNV.seqz.W, CNV.seqz)

# CNV
CNVGroupInfo <- get_groups(Sig.CNV, method = "consensus", match_consensus = TRUE)
CNVInfo <- CNV@summary.per.sample
CNVExposureInfo <- get_sig_exposure(Sig.CNV)
CNVscores <- scoring(CNV)

# Processing mutation data ------------------------------------------------

load(file = "output/PRAD_TCGA_plus_dbGap_Maf.RData")
load(file = "output/Sig.PRAD_TCGA_plus_dbGap_Maf.RData")

Maf_samp = data.table::data.table(
  Tumor_Sample_Barcode = unique(rbind(Maf@data, Maf@maf.silent)$Tumor_Sample_Barcode)
)

TMBInfo <- Maf@variant.type.summary
TMBInfo$n_INDEL = TMBInfo$DEL + TMBInfo$INS
TMBInfo = TMBInfo[, .(Tumor_Sample_Barcode, n_INDEL, SNP)]
colnames(TMBInfo)[3] = "n_SNV"
TMBInfo = merge(Maf_samp, TMBInfo, by = "Tumor_Sample_Barcode", all = TRUE)
# Fill NAs
TMBInfo = TMBInfo %>%
  dtplyr::lazy_dt() %>%
  dplyr::mutate_at(vars(dplyr::starts_with("n_")), ~ifelse(is.na(.), 0, .)) %>%
  data.table::as.data.table()

## Driver info
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

DriverDF = merge(Maf_samp, DriverDF, by.x = "Tumor_Sample_Barcode", by.y = "sample", all = TRUE)
# Fill NAs
DriverDF$n_driver = ifelse(is.na(DriverDF$n_driver), 0, DriverDF$n_driver)

save(DriverInfo, driver_genes, DriverDF, file = "output/PRAD_driver_info.RData")

##
TitvInfo <- titv(maf = Maf, plot = FALSE, useSyn = TRUE)$TiTv.fractions
MathInfo <- inferHeterogeneity(Maf, TitvInfo$Tumor_Sample_Barcode, useSyn = TRUE)
MathDF <- MathInfo$clusterData[, list(MATH = mean(MATH, na.rm = TRUE)), by = Tumor_Sample_Barcode]
ClusterDF <- MathInfo$clusterMeans[, list(cluster = ifelse(cluster=="outlier", 0L, as.integer(cluster)), Tumor_Sample_Barcode)][, list(cluster = max(cluster, na.rm = TRUE)), by = Tumor_Sample_Barcode]

# Fill NAs
MathDF = merge(Maf_samp, MathDF, by = "Tumor_Sample_Barcode", all = TRUE)
ClusterDF = merge(Maf_samp, ClusterDF, by = "Tumor_Sample_Barcode", all = TRUE)

fill_na = function(x, fill_value = 0) {
  x = ifelse(is.na(x), fill_value, x)
  x
}

summary(MathDF$MATH)
summary(ClusterDF$cluster)

MathDF$MATH = fill_na(MathDF$MATH)
ClusterDF$cluster = fill_na(ClusterDF$cluster, fill_value = 1) # Minimal should be 1 cluster

save(TitvInfo, MathInfo, MathDF, ClusterDF, file = "output/PRAD_heter_info.RData")

SNVGroupInfo <- get_groups(Sig.SNV, method = "consensus", match_consensus = TRUE)
SNVExposureInfo <- get_sig_exposure(Sig.SNV)

# Processing gene and pathway mutation ------------------------------------

## Extract gene variant info
gene_list1 = c("BRCA1", "BRCA2", "TP53", "SPOP", "FOXA1", "RB1", "PTEN", "AR",
               "APC", "CDK12", "PIK3CA", "PIK3CB", "ATM", "IDH1", "MED12")
# CTNNB1 check only oncogenic mutation which located in 32-45 of Protein

# AR: amp
# CDKN1B, NKX3-1, PIK3R1, PTEN, RB1, APC: del
gene_list2 = c("AR", "CDKN1B", "NKX3-1", "PIK3R1", "PTEN", "RB1", "APC")

tb_mutation = Maf@data[, .(Tumor_Sample_Barcode, Hugo_Symbol, Chromosome, Start_Position, End_Position,
                           Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2,
                           Protein_position, Consequence, SIFT, PolyPhen, CLIN_SIG)]
table(tb_mutation$Consequence)
# The amino acid substitution is predicted damaging is the score is <= 0.05, and tolerated if the score is > 0.05
table(tb_mutation$SIFT)
# > 0.15, if damaging
table(tb_mutation$PolyPhen)
#
table(tb_mutation$CLIN_SIG)

## Assume a mutation is functional based on SIFT/PolyPhen/CLIN_SIG
tb_mutation$is_pathogenic = grepl("deleterious", tb_mutation$SIFT) | grepl("damaging", tb_mutation$PolyPhen) | grepl("pathogenic", tb_mutation$CLIN_SIG)


# get_mutation_status <- function(x, genes) {
#   as_tibble(
#     setNames(
#       lapply(genes, function(gene) any(grepl(paste0("^", gene, "$"), x))),
#       genes
#     )
#   )
# }
# summary_mutation <- Maf@data %>%
#   group_by(Tumor_Sample_Barcode) %>%
#   do(get_mutation_status(.$Hugo_Symbol, driver_genes))

get_mutation_status2 <- function(x, y, genes) {
  as_tibble(
    setNames(
      lapply(genes, function(gene) any(grepl(paste0("^", gene, "$"), x) & y) ),
      genes
    )
  )
}

summary_mutation = tb_mutation %>%
  as_tibble() %>%
  group_by(Tumor_Sample_Barcode) %>%
  do(get_mutation_status2(.$Hugo_Symbol, .$is_pathogenic, gene_list1)) %>%
  data.table::as.data.table()

tb_mutation %>%
  dplyr::filter(Hugo_Symbol == "CTNNB1") %>%
  mutate(pos = as.integer(sub("([0-9]*-)?([0-9]+)/[0-9]+.", "\\2", Protein_position))) %>%
  select(Protein_position, pos)

mut_CTNNB1 = tb_mutation %>%
  as_tibble() %>%
  dplyr::filter(Hugo_Symbol == "CTNNB1") %>%
  mutate(pos = as.integer(sub("([0-9]*-)?([0-9]+)/[0-9]+.", "\\2", Protein_position))) %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(
    CTNNB1 = any(is_pathogenic | (pos >= 32 & pos <= 45 & !grepl("frameshift_variant", Consequence)))
  ) %>%
  data.table::as.data.table()

summary_mutation = merge(summary_mutation, mut_CTNNB1, by = "Tumor_Sample_Barcode", all = TRUE)

# Process CNV for genes
# GISTIC2 data
all.lesions <- "data/all_lesions.conf_99_v2.txt"
amp.genes <- "data/amp_genes.conf_99.txt"
del.genes <- "data/del_genes.conf_99.txt"
scores.gis <- "data/scores.gistic"

gistic = readGistic(gisticAllLesionsFile = all.lesions,
                    gisticAmpGenesFile = amp.genes,
                    gisticDelGenesFile = del.genes,
                    gisticScoresFile = scores.gis)

get_cnv_status <- function(x, y, genes, types) {
  as_tibble(
    setNames(
      purrr::map2(genes, types, function(gene, type) {
        any(grepl(paste0("^", gene, "$"), x) & y == type)
      }),
      genes
      )
    )
}


gene_list2
type_list = c("Amp", rep("Del", 6))

tb_cnv = gistic@data[Hugo_Symbol %in% gene_list2] %>%
  as_tibble() %>%
  group_by(Tumor_Sample_Barcode) %>%
  do(get_cnv_status(.$Hugo_Symbol, .$Variant_Classification,
                    gene_list2, type_list)) %>%
  ungroup() %>%
  mutate(Tumor_Sample_Barcode = as.character(Tumor_Sample_Barcode)) %>%
  data.table::as.data.table()

# Fill NAs and merge mutation for AR, PTEN, APC, RB1
summary_mutation = merge(Maf_samp, summary_mutation, by = "Tumor_Sample_Barcode", all = TRUE)
summary_mutation = merge(summary_mutation, tb_cnv, by = "Tumor_Sample_Barcode", all.x = TRUE)
summary_mutation = summary_mutation %>%
  as_tibble() %>%
  mutate_at(vars(-Tumor_Sample_Barcode), fill_na, fill_value = FALSE) %>%
  mutate(
    RB1 = RB1.x | RB1.y,
    PTEN = PTEN.x | PTEN.y,
    AR = AR.x | AR.y,
    APC = APC.x | APC.y
  ) %>%
  select(-c(paste0("RB1", c(".x", ".y")),
            paste0("PTEN", c(".x", ".y")),
            paste0("AR", c(".x", ".y")),
            paste0("APC", c(".x", ".y")))) %>%
  data.table::as.data.table()

# Only functional variants included here
## get mutation status of pathways

# get_pathway_mut_status <- function(x, paths) {
#   as_tibble(
#     setNames(
#       lapply(paths, function(path) any(grepl(paste0("^", path, "$"), x))),
#       paths
#     )
#   )
# }

get_pathway_mut_status2 <- function(x, y, paths) {
  as_tibble(
    setNames(
      lapply(paths, function(path) any(grepl(paste0("^", path, "$"), x) & y) ),
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

# AR pathway from GO MF: http://www.informatics.jax.org/go/term/GO:0030521
AR = readr::read_tsv("data/GO_MF_AR.txt")
AR = toupper(unique(AR$Symbol))

AR <- data.table::data.table(
  Pathway = "AR_pathway",
  Gene = AR
)

pathdb <- rbind(pathdb, HR, AR, fill = TRUE)
readr::write_csv(pathdb, path = "output/pathway_list.csv")

path_list <- unique(pathdb$Pathway)
mut_data <- dplyr::left_join(tb_mutation %>% as_tibble,
                             pathdb %>% as_tibble,
                             by = c("Hugo_Symbol" = "Gene"))

summary_pathway <- mut_data %>%
  dplyr::group_by(Tumor_Sample_Barcode) %>%
  dplyr::do(get_pathway_mut_status2(.$Pathway, .$is_pathogenic, path_list))

# Set pathway postfix
colnames(summary_pathway) = ifelse(
  endsWith(colnames(summary_pathway), "_pathway") | colnames(summary_pathway) == "Tumor_Sample_Barcode",
  colnames(summary_pathway),
  paste0(colnames(summary_pathway), "_pathway")
)

summary_pathway = merge(Maf_samp, summary_pathway, by = "Tumor_Sample_Barcode", all = TRUE)
summary_pathway = summary_pathway %>%
  as_tibble() %>%
  mutate_at(vars(-Tumor_Sample_Barcode), fill_na, fill_value = FALSE) %>%
  data.table::as.data.table()


save(summary_mutation, summary_pathway, file = "output/PRAD_gene_and_pathway_mutation.RData")

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

colnames(CNVExposureInfo) <- c("sample", paste0("CNV-", colnames(CNVExposureInfo)[-1]))
colnames(TMBInfo)[1] <- "sample"
TMBInfo[, sample := as.character(sample)][, sample := ifelse(startsWith(sample, "TCGA"),
                                                             substr(sample, 1, 15),
                                                             sample
)]

colnames(DriverDF)[1] = "sample"
DriverDF[, sample := as.character(sample)][, sample := ifelse(startsWith(sample, "TCGA"),
                                                              substr(sample, 1, 15),
                                                              sample
)]
colnames(SNVGroupInfo) <- c("sample", "snv_group", "snv_weight", "snv_enrich_sig")
SNVGroupInfo[, sample := as.character(sample)][, sample := ifelse(startsWith(sample, "TCGA"),
                                                                  substr(sample, 1, 15),
                                                                  sample
)]
colnames(SNVExposureInfo) <- c("sample", paste0("SBS-", colnames(SNVExposureInfo)[-1]))
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
colnames(ClusterDF) <- c("sample", "n_mutation_cluster")
ClusterDF[, sample := as.character(sample)][, sample := ifelse(startsWith(sample, "TCGA"),
                                                               substr(sample, 1, 15),
                                                               sample
)]

# data.table::setDT(summary_mutation)
# data.table::setDT(summary_pathway)

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


colnames(CNVInfo)[2:4] = c("n_cnv", "n_amp", "n_del")

MergeInfo <- Info[which(!is.na(Info$tumor_Run)), ] %>%
  left_join(CNVInfo[, .(sample, n_cnv, n_amp, n_del)], by = c("CNV_ID" = "sample")) %>%
  left_join(CNVscores, by = c("CNV_ID" = "sample")) %>%
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


saveRDS(MergeInfo, file = "output/PRAD_Merge_Info_CNV_from_sequenza_update.RDS")
