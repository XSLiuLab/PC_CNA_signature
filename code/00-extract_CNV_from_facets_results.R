# NOTE: the raw data are big, thus not maintained in repository
library(dplyr)
library(stringr)

# Obtaining absolute copy number from FACETS ------------------------------

extract_facets_cnv("raw_data/facets_results/wgs/", target_path = "data/CNV_from_TCGA_WGS.tsv")
extract_facets_cnv("raw_data/facets_results/wes/", target_path = "data/CNV_from_TCGA_WES.tsv")
# //TODO: currently, samples from dbGap are not processed by this function
