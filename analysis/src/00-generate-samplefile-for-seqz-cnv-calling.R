# FUN: generate sample file for sequenza copy number calling
library(tidyverse)

load("dbGap/mapping_df.RData")
TCGA_Sample = data.table::fread("data/TCGA_sample_id.txt", header = FALSE)

get_path = function(study, name) {
    # Proj 21926 contains
    Proj21926 = paste0("phs", c("000915", "001141", "000554"))
    Proj16533 = paste0("phs", c("000447", "000909"))
    Proj_TCGA = "TCGA"
    ifelse(
        study %in% Proj_TCGA,
        paste0("/public/home/liuxs/biodata/gdc/links/TCGA_PRAD/",
                name, "*.bam"),
        ifelse(
            study %in% Proj16533,
            paste0("/public/home/liuxs/ncbi/dbGaP-16533/dnaseq/BQSR/bqsrbam/", name, ".sorted.marked.BQSR.bam"),
            ifelse(
                study %in% Proj21926, 
                paste0("/public/home/liuxs/ncbi/dbGaP-21926/dnaseq/BQSR/bqsrbam/",name, ".sorted.marked.BQSR.bam"),
                NA_character_)
        )
    )
}

All_data = dplyr::bind_rows(
    mapping_df %>% 
        dplyr::rename(Study = gap_accession) %>%
        dplyr::mutate(sample = paste(subject_id, tumor_Run, sep="-")) %>%
        dplyr::select(Study, sample, tumor_Run, normal_Run),
    TCGA_Sample %>%
        rename(tumor_Run=V1, normal_Run=V2) %>%
        dplyr::mutate(Study = "TCGA",
                     sample = sub(".*(TCGA-[^\\.]+).*", "\\1", tumor_Run)) %>% 
        dplyr::mutate(sample = substr(sample, 1, 15))

)

WES_Info = All_data %>% 
    dplyr::mutate(
        tumor_Run = get_path(Study, tumor_Run),
        normal_Run = get_path(Study, normal_Run)
    ) %>% 
    dplyr::select(-Study)

readr::write_csv(WES_Info, path="cnv_calling/sequenza/pipeline/WES_samples.csv", col_names=FALSE)