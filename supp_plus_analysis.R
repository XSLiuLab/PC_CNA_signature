library(tidyverse)

# 01 检查是否存在一个患者的原发和转移样本 ---------------------------------------------------

df.seqz = readRDS(file = "output/PRAD_Merge_Info_CNV_from_sequenza_update.RDS")
df_dup <- df.seqz[duplicated(df.seqz$subject_id), ]
df_dup %>% filter(sample_type == "Primary") %>% pull(subject_id)

View(df_dup %>% filter(subject_id %in% c("WCMC7520","WCMC9")))

df_dup %>% filter(subject_id == "WCMC7520") %>%
  select(subject_id, tumor_body_site, `CN-Sig1`, `CN-Sig2`, `CN-Sig3`, `CN-Sig4`, `CN-Sig5`) %>%
  mutate(across(where(is.numeric), ~round(., 2)))


# 02 是否有与 ecDNA 研究交叉的样本 ---------------------------------------------------
library(readxl)

tbl1 <- read_excel("data/ICGC-ecDNA.xlsx", sheet = 1)
ecDNA <- tbl1 %>%
  count(sample_barcode, amplicon_classification)

common_names <- intersect(ecDNA$sample_barcode, df.seqz$subject_id)
View(df.seqz %>% filter(subject_id %in% common_names))

df_ecDNA <- df.seqz %>% filter(subject_id %in% common_names)
tbl_ecDNA <- ecDNA %>% filter(sample_barcode %in% common_names)

df_ecDNA <- left_join(tbl_ecDNA,
          df_ecDNA %>% select(subject_id, `CN-Sig1`:`CN-Sig5`),
          by = c("sample_barcode" = "subject_id"))

df_ecDNA <- df_ecDNA %>%
  mutate(ac = glue::glue("{amplicon_classification}({n})")) %>%
  group_by(sample_barcode) %>%
  summarise(
    ac = paste(ac, collapse = "+")
  ) %>%
  left_join(df_ecDNA %>% select(sample_barcode, `CN-Sig1`:`CN-Sig5`) %>% unique()) %>%
  arrange(desc(ac))

library(ggpubr)

p <- ggbarplot(df_ecDNA %>% rename(amplicon_type = ac), x = "sample_barcode", y = c("CN-Sig1"),
          sort.by.groups = FALSE, fill = "amplicon_type", label = df_ecDNA$ac, xlab = FALSE) + rotate_x_text()
ggsave("output/amplicon_and_CNS1.png", plot = p, width = 9, height = 6)

p <- ggbarplot(df_ecDNA %>% rename(amplicon_type = ac), x = "sample_barcode", y = c("CN-Sig2"),
               sort.by.groups = FALSE, fill = "amplicon_type", label = df_ecDNA$ac, xlab = FALSE) + rotate_x_text()
ggsave("output/amplicon_and_CNS2.png", plot = p, width = 9, height = 6)


# 03 Match FACETS and Sequenza --------------------------------------------

library(sigminer)
library(data.table)
load("output/CNV.facets.RData")
load("output/CNV.seqz.RData")

df.seqz = readRDS(file = "output/PRAD_Merge_Info_CNV_from_sequenza_update.RDS") %>%
  select(CNV_ID, tumor_Run, normal_Run) %>%
  rename(id_seqz = CNV_ID)
df.facets = readRDS(file = "output/PRAD_Merge_Info_CNV_from_facets.RData") %>%
  select(CNV_ID, tumor_Run, normal_Run) %>%
  rename(id_facets = CNV_ID) %>% na.omit()

id_map <- full_join(df.seqz, df.facets)
id_map2 <- id_map$id_seqz
names(id_map2) <- id_map$id_facets

df = data.table::fread("data/CNV_from_dbGAP_PLUS_TCGA_WES_CVAL150.tsv")
df$sample <- id_map2[df$sample]
df = df[sample != "10-SRR4043970"]
df = df %>%
  mutate(Chromosome = paste0("chr", Chromosome),
         Chromosome = case_when(
           Chromosome == "chr23" ~ "chrX",
           Chromosome == "chr24" ~ "chrY",
           TRUE ~ Chromosome
         ))

df2 = data.table::fread("data/CNV_from_sequenza.tsv")
all(df$sample %in% df2$sample)

#========
df$facets_idx = 1:nrow(df)
df2$seqz_idx = 1:nrow(df2)

all_samps = unique(df2$sample)

result = data.table()
for (i in all_samps) {
  message("Handling sample ", i)
  z1 = df[sample == i]
  z2 = df2[sample == i]
  setkey(z2, Chromosome, Start.bp, End.bp)
  z = foverlaps(z1, z2)
  result = rbind(result, z)
}

result$i.modal_cn = ifelse(result$i.modal_cn > 20, 20, result$i.modal_cn)

result2 = na.omit(result)
result_sim <- result2 %>%
  group_by(i.sample) %>%
  summarise(sim = cosine(modal_cn, i.modal_cn))

pdf("figures_new/FACETS_Sequenza_similarity.pdf", width = 5, height = 4)
plot(density(result_sim$sim), main = "Copy number profile similarity\nfrom two methods")
dev.off()

pdf("figures_new/FACETS_Sequenza_similarity_hist.pdf", width = 5, height = 4)
hist(result_sim$sim,
     main = "Copy number profile similarity\nfrom two methods", xlab = NA)
dev.off()

# 找到和对比删除长度
del1 = df[modal_cn == 0 & !Chromosome %in% c("chrX", "chrY")]
del2 = df2[modal_cn == 0 & !Chromosome %in% c("chrX", "chrY")]
del1$Len = (del1$End.bp - del1$Start.bp + 1L) / 1e6
del2$Len = (del2$End.bp - del2$Start.bp + 1L) / 1e6

# del1 %>%
#   group_by(sample) %>%
#   summarise(count = sum(Len > 10)) %>%
#   arrange(desc(count)) -> del1_s
#
# del2 %>%
#   group_by(sample) %>%
#   summarise(count = sum(Len > 10)) %>%
#   arrange(desc(count)) -> del2_s
#
# full_join(del1_s, del2_s, by = "sample") %>%
#   mutate(df = count.x - count.y) %>%
#   arrange(count.y, desc(df))

#==
del1 %>%
  group_by(sample) %>%
  summarise(Len = sum(Len)) %>%
  arrange(desc(Len)) -> del1_s

del2 %>%
  group_by(sample) %>%
  summarise(Len = sum(Len)) %>%
  arrange(desc(Len)) -> del2_s

full_join(del1_s, del2_s, by = "sample") %>%
  mutate(df = Len.x - Len.y) %>%
  left_join(result_sim, by = c("sample"="i.sample")) %>%
  arrange(desc(df), Len.y, sim)

p1 <- show_cn_profile(
  df %>% as_tibble() %>% setnames(c("chromosome", "start", "end", "segVal", "sample", "idx")),
  samples = c("09-5446-SRR473087", "1115081-SRR8304753", "TCGA-HC-7750-01", "5115056-SRR8311885"), show_title = TRUE)

p2 <- show_cn_profile(
  df2 %>% as_tibble() %>%  setnames(c("chromosome", "start", "end", "segVal", "sample", "idx")),
  samples = c("09-5446-SRR473087", "1115081-SRR8304753", "TCGA-HC-7750-01", "5115056-SRR8311885"), show_title = TRUE)

ggsave(filename = "figures_new/facets_del_sample_profile.pdf", p1,
       height = 6, width = 12)
ggsave(filename = "figures_new/seqz_facets_del_sample_profile.pdf", p2,
       height = 6, width = 12)



###
library(sigminer)
library(tidyverse)
library(NMF)

tcga_prad = readRDS("tcga_prad.rds") %>%
  mutate(sample = substr(sample_submitter_id, 1, 15),
         Start.bp = Start,
         End.bp = End,
         modal_cn = Copy_Number) %>%
  select(Chromosome, Start.bp, End.bp, modal_cn, sample)

# Set this per R session
options(sigminer.sex = "male", sigminer.copynumber.max = 20L)

CNV.tcga_prad <- read_copynumber(
  tcga_prad,
  genome_build = "hg38",
  complement = FALSE, verbose = TRUE
)

tally.W <- sig_tally(CNV.tcga_prad, method = "W", cores = 6, feature_setting = CN.features)

# EST <- sig_estimate(tally.W$nmf_matrix,
#                     range = 2:6, nrun = 50,
#                     cores = 1,
#                     save_plots = FALSE,
#                     verbose = TRUE
# )
bp_prad <- bp_extract_signatures(tally.W$nmf_matrix, range = 2:6, n_bootstrap = 5, n_nmf_run = 10, cores = 4, pynmf = TRUE)
Sig.prad <- sig_extract(tally.W$nmf_matrix, n_sig = 4, nrun = 50, cores = 1)

load("output/Sig.CNV.seqz.W.RData")
sim <- get_sig_similarity(Sig.prad, Sig.CNV.seqz.W)

pheatmap::pheatmap(sim$similarity, cluster_rows = F, cluster_cols = F, display_numbers = TRUE, filename = "figures_new/tcga_prad_denovo.pdf", height = 4, width = 4.5)

# TCGA 拷贝数图谱相关性
all_samps = unique(tcga_prad$sample)
df2 = data.table::fread("data/CNV_from_sequenza.tsv")

result = data.table()
for (i in all_samps) {
  message("Handling sample ", i)
  z1 = tcga_prad[sample == i]
  z2 = df2[sample == i]
  setkey(z2, Chromosome, Start.bp, End.bp)
  z = foverlaps(z1, z2)
  result = rbind(result, z)
}

result$i.modal_cn = ifelse(result$i.modal_cn > 20, 20, result$i.modal_cn)

result2 = na.omit(result)
result_sim <- result2 %>%
  group_by(i.sample) %>%
  summarise(sim = cosine(modal_cn, i.modal_cn))

pdf("figures_new/TCGA_PRAD_SNP_and_sequenza_WES_CNP_similarity.pdf", width = 5, height = 4)
plot(density(result_sim$sim), main = "Copy number profile similarity")
dev.off()
