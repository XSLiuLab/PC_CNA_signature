# phs000554 should also be used in survival analysis

library(UCSCXenaTools)  # >= 1.2.10
library(tidyverse)
library(dtplyr)
library(data.table)

x_query = XenaData %>%
  XenaGenerate(subset = XenaHostNames == "tcgaHub") %>%
  XenaFilter(filterDatasets = "PRAD") %>%
  XenaFilter(filterDatasets = "clinical|survival") %>%
  XenaQuery()
x_download = x_query %>%
  XenaDownload(destdir = "data/UCSCXena", trans_slash = TRUE)

TCGA_cli = data.table::fread(cli_download$destfiles)

# Curated clinical data from the Pan-cancer Atlas paper titled "An Integrated TCGA Pan-Cancer Clinical Data Resource (TCGA-CDR) to drive high quality survival outcome analytics". The paper highlights four types of carefully curated survival endpoints, and recommends the use of the endpoints of OS, PFI, DFI, and DSS for each TCGA cancer type.
#
# OS: overall survial
# PFI: progression-free interval
# DSS: disease-specific survival
# DFI: disease-free interval

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
