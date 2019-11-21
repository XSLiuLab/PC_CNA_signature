# Extract samples with deep deletion feature 
# to do GISTIC analysis
# for finding recurrent regions
library(sigminer)
library(dtplyr)

seg = readr::read_tsv("data/PRAD_CNA_hg38.seg")
CNV.Sig = readRDS("output/NMF_copynumber_signature.prob.rds")


rel_expo = get_sig_exposure(CNV.Sig, type = "relative")
rel_expo = lazy_dt(rel_expo)
samps = rel_expo %>%
    dplyr::filter(Sig3 > 0.2) %>%
    dplyr::pull(sample)

data = seg %>% 
    dplyr::filter(ID %in% samps)

readr::write_tsv(data, "data/PRAD_CNA_hg38_for_sig3.seg")