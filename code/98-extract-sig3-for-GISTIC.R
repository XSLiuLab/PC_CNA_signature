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

# The GISTIC2 results for samples above is chaos
# So we try to focus on Sig3 dominant samples
groups = get_groups(CNV.Sig)
samps2 = groups %>% 
    lazy_dt() %>% 
    dplyr::filter(enrich_sig == "Sig3") %>% 
    dplyr::pull(sample)

data2 = seg %>% 
    dplyr::filter(ID %in% samps2)

readr::write_tsv(data2, "data/PRAD_CNA_hg38_for_sig3_dominant.seg")