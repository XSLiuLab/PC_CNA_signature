library(tidyverse)
library(openxlsx)

pathway = readr::read_csv("output/pathway_list.csv")
pathway = pathway %>%
  dplyr::select(-OG_TSG) %>%
  dplyr::mutate(Pathway = ifelse(endsWith(Pathway, "_pathway"),
                                 sub("_pathway", "", Pathway),
                                 Pathway))

BAM_summary = readr::read_csv("output/PRAD_PROJ_BAM_SUMMARY.csv", col_names = FALSE)
colnames(BAM_summary) = c("Path", "Count", "Reads")

BAM_summary = BAM_summary %>%
  dplyr::mutate(Average = Reads / Count) %>%
  dplyr::mutate(Name = basename(Path),
                Name = sub("[_-]depth", "", Name),
                Name = ifelse(grepl("TCGA", Name),
                              substr(Name, 1+5, 15+5),
                              Name)) %>%
  dplyr::select(-Path) %>%
  dplyr::distinct(Name, .keep_all = TRUE)

df = readRDS("output/df.seqz.RDS")

df_merge = df %>%
  dplyr::select(tumor_Run, normal_Run) %>%
  tidyr::pivot_longer(cols = c("tumor_Run", "normal_Run"), names_to = "Type", values_to = "ID") %>%
  dplyr::left_join(BAM_summary, by = c("ID"="Name"))

df_bam = df %>%
  dplyr::select(tumor_Run, normal_Run) %>%
  dplyr::left_join(
    dplyr::filter(df_merge, Type == "tumor_Run") %>% dplyr::select(-Type) %>% unique(),
    by = c("tumor_Run"="ID")
  ) %>%
  setNames(c("tumor_Run", "normal_Run", "Count (depth > 25)", "Reads (depth > 25)", "Average (depth > 25)")) %>%
  dplyr::left_join(
    dplyr::filter(df_merge, Type == "normal_Run") %>% dplyr::select(-Type) %>% unique(),
    by = c("normal_Run"="ID")
  ) %>%
  unique()

df_bam = df_bam[, c(1,3:5, 2, 6:8)]
colnames(df_bam) = c("Tumor ID", "Tumor counts", "Tumor reads", "Tumor average depth",
                     "Normal ID", "Normal counts", "Tumor reads", "Normal average depth")

summary(df_bam$`Tumor average depth`)
summary(df_bam$`Normal average depth`)

write.xlsx(
  list(df_bam, pathway), file = "output/BAM_and_gene_list.xlsx"
)


load("output/CNV.seqz.RData")
write.xlsx(
  CNV.seqz@data,
  file = "output/CNV_seginfo.xlsx"
)
