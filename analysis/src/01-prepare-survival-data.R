# Curated clinical data from the Pan-cancer Atlas paper titled "An Integrated TCGA Pan-Cancer Clinical Data Resource (TCGA-CDR) to drive high quality survival outcome analytics". The paper highlights four types of carefully curated survival endpoints, and recommends the use of the endpoints of OS, PFI, DFI, and DSS for each TCGA cancer type.
#
# OS: overall survial
# PFI: progression-free interval
# DSS: disease-specific survival
# DFI: disease-free interval
#
# phs000554 should also be used in survival analysis

library(UCSCXenaTools) # >= 1.2.10
library(tidyverse)
library(dtplyr)
library(data.table)


# TCGA --------------------------------------------------------------------

x_query <- XenaData %>%
  XenaGenerate(subset = XenaHostNames == "tcgaHub") %>%
  XenaFilter(filterDatasets = "PRAD") %>%
  XenaFilter(filterDatasets = "clinical|survival") %>%
  XenaQuery()
x_download <- x_query %>%
  XenaDownload(destdir = "data/UCSCXena", trans_slash = TRUE)


tcga_prad_cli <- XenaPrepare(x_download)
tcga_prad_surv <- tcga_prad_cli$survival__PRAD_survival.txt.gz %>%
  dplyr::select(-`_PATIENT`, -Redaction) %>%
  dplyr::mutate(Study = "TCGA")




# phs000554 ---------------------------------------------------------------

PRAD_cli <- readRDS(file = "data/PRAD_CLINICAL.rds")
phs000554 <- readRDS("data/Tidy_Clinical/phs000554.rds")

phs000554_surv <- phs000554 %>%
  dplyr::rename(OS.time = `Survival from diagnosis (mo)`) %>%
  dplyr::filter(!is.na(OS.time)) %>%
  dplyr::select(SampleName, OS.time) %>%
  dplyr::rename(sample = SampleName) %>%
  dplyr::mutate(sample = ifelse(sample == "WA43-44 (bladder)", "WA43", sample)) %>%
  dplyr::mutate(OS = 1L) %>% # No dead status reported, so assume they are dead due to their metastatic status
  dplyr::mutate(
    Study = "phs000554",
    OS.time = OS.time * 30L
  )



# Merge -------------------------------------------------------------------

prad_surv <- dplyr::bind_rows(tcga_prad_surv, phs000554_surv)
saveRDS(prad_surv, file = "data/PRAD_Survival.rds")


# PRAD_cli = TCGA_cli %>%
#   lazy_dt() %>%
#   dplyr::filter(sample %in% TCGA_IDs) %>%
#   dplyr::select(sample, DFI, DFI.time, DSS, DSS.time, OS, OS.time, PFI, PFI.time) %>%
#   dplyr::mutate(sample = substr(sample, 1, 12)) %>%
#   as.data.table()
#
# summary(PRAD_cli)
#
# load(file = "data/PRAD_Merge_Info.RData")
#
# TCGA_Df = MergeInfo %>%
#   dplyr::filter(Study == "TCGA") %>%
#   dplyr::inner_join(PRAD_cli, by = c("PatientID"="sample"))
#
#
# library(ezcox)
# ezcox(TCGA_Df, covariates = c(paste0("CNV_Sig", 1:6), paste0("SNV_Sig", 1:3)), time = "OS.time", status = "OS")
# ezcox(TCGA_Df, covariates = c(paste0("CNV_Sig", 1:6), paste0("SNV_Sig", 1:3)), time = "PFI.time", status = "PFI")
# ezcox(TCGA_Df, covariates = c(paste0("CNV_Sig", 1:6), paste0("SNV_Sig", 1:3)), time = "DSS.time", status = "DSS")
# ezcox(TCGA_Df, covariates = c(paste0("CNV_Sig", 1:6), paste0("SNV_Sig", 1:3)), time = "DFI.time", status = "DFI")
#
# zz = ezcox(TCGA_Df, covariates = c(paste0("CNV_Sig", 1:6), paste0("SNV_Sig", 1:3)), time = "OS.time", status = "OS", return_models = TRUE)
# zz$models$model
# ?forest_model()
