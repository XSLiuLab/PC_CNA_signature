# Extract samples with deep deletion feature
# to do GISTIC analysis
# for finding recurrent regions
library(sigminer)
library(dtplyr)

seg <- readr::read_tsv("data/PRAD_CNA_hg38.seg")
CNV.Sig <- readRDS("output/NMF_copynumber_signature.prob.rds")


rel_expo <- get_sig_exposure(CNV.Sig, type = "relative")
rel_expo <- lazy_dt(rel_expo)
samps <- rel_expo %>%
  dplyr::filter(Sig3 > 0.2) %>%
  dplyr::pull(sample)

data <- seg %>%
  dplyr::filter(ID %in% samps)

readr::write_tsv(data, "data/PRAD_CNA_hg38_for_sig3.seg")

# The GISTIC2 results for samples above is chaos
# So we try to focus on Sig3 dominant samples
groups <- get_groups(CNV.Sig)
samps2 <- groups %>%
  lazy_dt() %>%
  dplyr::filter(enrich_sig == "Sig3") %>%
  dplyr::pull(sample)

data2 <- seg %>%
  dplyr::filter(ID %in% samps2)

readr::write_tsv(data2, "data/PRAD_CNA_hg38_for_sig3_dominant.seg")

# Generate sample file for sequenza calling
Info <- readRDS("data/PRAD_CLINICAL.rds")
Del_Info <- Info %>%
  dplyr::filter(CNV_ID %in% samps2) %>%
  dplyr::select(CNV_ID, tumor_Run, normal_Run, Study)

get_path <- function(study, name) {
  # Proj 21926 contains
  Proj21926 <- paste0("phs", c("000915", "001141", "000554"))
  Proj16533 <- paste0("phs", c("000447", "000909"))
  Proj_TCGA <- "TCGA"
  ifelse(
    study %in% Proj_TCGA,
    paste0(
      "/public/home/liuxs/biodata/gdc/links/TCGA_PRAD/",
      name, "*.bam"
    ),
    ifelse(
      study %in% Proj16533,
      paste0("/public/home/liuxs/ncbi/dbGaP-16533/dnaseq/BQSR/bqsrbam/", name, ".sorted.marked.BQSR.bam"),
      ifelse(
        study %in% Proj21926,
        paste0("/public/home/liuxs/ncbi/dbGaP-21926/dnaseq/BQSR/bqsrbam/", name, ".sorted.marked.BQSR.bam"),
        NA_character_
      )
    )
  )
}

TCGA_Sample <- data.table::fread("data/TCGA_sample_id.txt", header = FALSE)

Del_Info$tumor_Run <- sapply(Del_Info$tumor_Run, function(x) {
  if (startsWith(x, "TCGA")) {
    grep(x, TCGA_Sample$V1, value = TRUE)
  } else {
    x
  }
}) %>% as.character()

Del_Info$normal_Run <- sapply(Del_Info$normal_Run, function(x) {
  if (startsWith(x, "TCGA")) {
    grep(x, TCGA_Sample$V2, value = TRUE)
  } else {
    x
  }
}) %>% as.character()

Del_Info <- Del_Info %>%
  dplyr::mutate(
    tumor_Run = get_path(Study, tumor_Run),
    normal_Run = get_path(Study, normal_Run)
  ) %>%
  dplyr::select(-Study)

readr::write_csv(Del_Info, path = "cnv_calling/sequenza/pipeline/deep_del_samples.csv", col_names = FALSE)
