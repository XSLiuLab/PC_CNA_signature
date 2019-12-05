# NOTE: the raw data are big, thus not maintained in repository
library(dplyr)
library(stringr)
#source("code/99-functions.R")

# Obtaining absolute copy number from Sequenza ------------------------------
extract_seqz_cnv = function(target_dir, target_path) {
    SAMPLE = dir(target_dir, pattern = "_segments.txt") %>%
        stringr::str_remove("_segments.txt")
    res = purrr::map2_df(file.path(target_dir, 
        paste(SAMPLE, 
        "_segments.txt", 
        sep = "")), SAMPLE, function(x, y) {
            message("Processing ", y)
            df = readr::read_tsv(x)
            df = df %>% 
                dplyr::select(chromosome, start.pos, end.pos, CNt) %>%
                dplyr::mutate(sample = y)
            colnames(df) = c("Chromosome","Start.bp","End.bp","modal_cn","sample")
            df
        })
    readr::write_tsv(res, path=target_path)
}

extract_seqz_cnv("/home/wsx/projects/prad_signature/cnv_calling/sequenza/seqz_wes_result",
                   target_path = "data/CNV_from_sequenza.tsv")

# Test
seqz_cnv = readr::read_tsv("data/CNV_from_sequenza.tsv")
test = sigminer::read_copynumber(seqz_cnv,
    genome_build = "hg38",
    complement = FALSE, verbose = TRUE)

# Obtaining purity and ploidy ---------------------------------------------

extract_seqz_purity_and_ploidy = function(target_dir, target_path) {
    SAMPLE = dir(target_dir, pattern = "_alternative_solutions.txt") %>%
        stringr::str_remove("_alternative_solutions.txt")
    res = purrr::map2_df(file.path(target_dir, 
        paste(SAMPLE, 
        "_alternative_solutions.txt", 
        sep = "")), SAMPLE, function(x, y) {
            message("Processing ", y)
            df = readr::read_tsv(x)
            df = df %>% 
                dplyr::select(cellularity,ploidy) %>%
                dplyr::mutate(sample = y) %>%
                dplyr::slice(1)  # Choose the best solution
            colnames(df) = c("purity", "ploidy", "sample")
            df
        })
    readr::write_tsv(res, path=target_path)
}

extract_seqz_purity_and_ploidy("/home/wsx/projects/prad_signature/cnv_calling/sequenza/seqz_wes_result",
                                 target_path = 'data/PRAD_Purity_and_Ploidy_Sequenza.tsv')

# Transform data to GISTIC input
# and store as .seg file
facets_to_GISTIC2("/public/data/facet/dbGAP_PLUS_TCGA_PRAD_WES_CVAL150/",
                  target_path = 'data/PRAD_CNA_hg38.seg',
                  rm_samps = "141-10")
