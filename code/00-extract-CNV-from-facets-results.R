# NOTE: the raw data are big, thus not maintained in repository
library(dplyr)
library(stringr)
source("code/99-functions.R")

# Obtaining absolute copy number from FACETS ------------------------------

extract_facets_cnv("/public/data/facet/dbGAP_PLUS_TCGA_PRAD_WES_CVAL150/",
                   target_path = "data/CNV_from_dbGAP_PLUS_TCGA_WES_CVAL150.tsv")

extract_facets_cnv("/public/data/facet_300/facetdata_300/",
                   target_path = "data/CNV_from_dbGAP_PLUS_TCGA_WES_CVAL300.tsv")

#extract_facets_cnv("raw_data/facets_results/wgs/", target_path = "data/CNV_from_TCGA_WGS.tsv")
#extract_facets_cnv("raw_data/facets_results/wes/", target_path = "data/CNV_from_TCGA_WES.tsv")


# Obtaining purity and ploidy ---------------------------------------------

extract_facets_purity_and_ploidy("/public/data/facet/dbGAP_PLUS_TCGA_PRAD_WES_CVAL150/",
                                 target_path = 'data/PRAD_Purity_and_Ploidy_CVAL150.tsv')

# Transform data to GISTIC input
# and store as .seg file
facets_to_GISTIC2("/public/data/facet/dbGAP_PLUS_TCGA_PRAD_WES_CVAL150/",
                  target_path = 'data/PRAD_CNA_hg38.seg')
