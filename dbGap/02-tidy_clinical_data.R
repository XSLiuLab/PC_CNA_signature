library(tidyverse)
library(readxl)
library(magrittr)

load("dbGap/dbGap_phenotype.RData")


# Clean datasets one by one -----------------------------------------------

# ====== phs000447 =======
phs000447 <- res_list$phs000447
phs000447_T1 <- read_excel(path = "data/TIdy_Clinical_Excel/phs000447.xlsx", sheet = 1)
phs000447_T2 <- read_excel(path = "data/TIdy_Clinical_Excel/phs000447.xlsx", sheet = 2)
phs000447_T3 <- read_excel(path = "data/TIdy_Clinical_Excel/phs000447.xlsx", sheet = 3)

phs000447_T1 <- phs000447_T1 %$%
  str_remove_all(PatientID, "-Tumor") %>%
  str_remove_all("-Normal") %>%
  unique()

phs000447_T3 <- phs000447_T3 %>%
  rename(
    SampleID = `Sample ID`
  ) %>%
  mutate(
    Sample = str_remove_all(Sample, "-Tumor") %>%
      str_remove_all("-Normal"),
    SampleID = gsub(SampleID, pattern = "(.+[0-9]+)[^0-9]*", replacement = "\\1") %>%
      gsub(pattern = "(.+[0-9]+[^_])_.+$", replacement = "\\1") %>%
      sub(pattern = "08-492T1", replacement = "08-492")
  ) %>%
  unique()

## Comment: T2 and T3 come from the same paper, so they have same samples
View(phs000447)

all(phs000447_T2$PatientID %in% phs000447_T3$Sample)
phs000447_T <- dplyr::full_join(phs000447_T3, phs000447_T2, by = c("Sample" = "PatientID"))
phs000447_T <- dplyr::full_join(
  dplyr::tibble(Sample = phs000447_T1),
  phs000447_T,
  by = "Sample"
)

# Why the samples have P or PR prefix???
# I cannot find the source from two original papers
phs000447$SUBJID[startsWith(phs000447$SUBJID, "STID")] %>% nchar()

phs000447_T$Predict_SampleID <- phs000447_T$Sample %>%
  str_remove("^PR-") %>%
  str_remove("^P") %>%
  {
    ifelse(str_detect(., "^[0-9]+$"), paste0("STID000000", .), .)
  }

# Sample IDs with STID prefix should all be 14
phs000447_T$Predict_SampleID[startsWith(phs000447_T$Predict_SampleID, "STID")] %>% nchar()
# Check if all these sample IDs exist in annotation data from dbGap
t1 <- phs000447 %>%
  dplyr::select(
    c(
      "dbGaP SubjID", "SUBJID", "Age", "Gender",
      "Primary Disease", "Cancer Stage", "PSA Pre-Operative",
      "Radiation Therapy"
    )
  ) %>%
  unique()
fill_value <- function(x) {
  if (all(is.na(x))) {
    return(NA_character_)
  } else {
    x_value <- x[!is.na(x)]
    if (length(unique(x_value)) > 1) {
      message("Inconsistent record detected, reset it to NA.")
      return(NA_character_)
    }
    return(x_value[1])
  }
}
t2 <- phs000447 %>%
  dplyr::select(c(
    "SUBJID", "Gleason Score",
    "TMPRSS2-ERG Fusion Status"
  )) %>%
  dplyr::group_by(SUBJID) %>%
  summarise(
    `Gleason Score` = fill_value(`Gleason Score`),
    `TMPRSS2-ERG Fusion Status` = fill_value(`TMPRSS2-ERG Fusion Status`)
  )

phs000447 <- dplyr::full_join(t1, t2, by = "SUBJID")
phs000447_T$Predict_SampleID[!phs000447_T$Predict_SampleID %in% phs000447$SUBJID]

phs000447 <- dplyr::full_join(phs000447_T %>% dplyr::select(Sample, Predict_SampleID),
  phs000447,
  by = c("Predict_SampleID" = "SUBJID")
)

# Check the rows cannot match
id_unmacth <- phs000447 %>%
  filter(is.na(`dbGaP SubjID`)) %>%
  filter(Predict_SampleID %in% phs000447_T$Predict_SampleID) %>%
  pull(Predict_SampleID)

phs000447_T %>%
  filter(Predict_SampleID %in% id_unmacth) %>%
  View()
# "03-022" has no extra information

# https://stackoverflow.com/questions/34096162/dplyr-mutate-replace-several-columns-on-a-subset-of-rows
mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}

t_supp <- phs000447_T %>%
  filter(Predict_SampleID %in% id_unmacth[-1])

phs000447 <- phs000447 %>%
  mutate_cond(Predict_SampleID %in% id_unmacth[-1],
    Age = t_supp$Age,
    `Primary Disease` = "Prostate Cancer",
    `Cancer Stage` = t_supp$`Pathological stage`,
    `PSA Pre-Operative` = t_supp$`Serum PSA at diagnosis (ng/mL)`,
    `Gleason Score` = t_supp$`Gleason Score`,
    `TMPRSS2-ERG Fusion Status` = t_supp$`ETS fusion detected by sequencing`
  ) %>%
  mutate(
    Age = as.integer(Age),
    `PSA Pre-Operative` = as.numeric(`PSA Pre-Operative`)
  )

rm(phs000447_T, phs000447_T1, phs000447_T2, phs000447_T3, t, t_supp, t1, t2, id_unmacth)
# dir.create("data/Tidy_Clinical")
saveRDS(phs000447, file = "data/Tidy_Clinical/phs000447.rds")

# ====== phs000554 =======
phs000554 <- res_list$phs000554
phs000554_T <- read_excel("data/TIdy_Clinical_Excel/phs000554.xlsx")

phs000554 <- phs000554 %>%
  filter(!grepl("matched", SAMPLE_ALIAS)) %>%
  select(-Sample_Name, -dbGaP_Sample_ID, -ANALYTE_TYPE, -IS_TUMOR)

phs000554 <- full_join(phs000554_T %>%
  select(-`Matched GE/aCGH`),
phs000554 %>%
  select(SUBJECT_ID, PRIMARY_METASTATIC_TUMOR, SAMPLE_ALIAS),
by = c("SampleName" = "SAMPLE_ALIAS")
) %>%
  mutate(
    Age = as.integer(Age),
    `Serum PSA` = as.numeric(`Serum PSA`),
    `Prior Treatment` = ifelse(`Prior Treatment` == "NA", NA, `Prior Treatment`)
  ) %>%
  mutate_at(vars(starts_with("Survival")), as.integer)

saveRDS(phs000554, file = "data/Tidy_Clinical/phs000554.rds")

# ====== phs000909 =======
phs000909 <- res_list$phs000909 # 这里面有部分全数字的编号，不知道CNV结果中有没有这些样本？
phs000909_T <- read_excel("data/TIdy_Clinical_Excel/phs000909.xlsx")
saveRDS(phs000909_T, file = "data/Tidy_Clinical/phs000909.rds")

# ====== phs000915 =======
phs000915 <- res_list$phs000915
phs000915_T <- read_excel("data/TIdy_Clinical_Excel/phs000915.xlsx")
all(phs000915$SUBJECT_ID %in% phs000915_T$`cBio_SU2C ID`)
which(phs000915$SUBJECT_ID %in% phs000915_T$`cBio_SU2C ID`) %>% length()

phs000915 <- full_join(
  phs000915_T %>%
    mutate(`cBio_SU2C ID` = as.character(`cBio_SU2C ID`)),
  phs000915 %>%
    select(SUBJECT_ID) %>%
    unique(),
  by = c("cBio_SU2C ID" = "SUBJECT_ID")
) %>%
  filter(!is.na(`cBio_SU2C ID`))
# 这个数据集文献只有 150 个样本，dbGap 结果显示为 146 个，
# 能够匹配的是 100 多点，
# 但 Nat.Gen 2018 文章却显示有 287 个
# 需要认真关注下
saveRDS(phs000915, file = "data/Tidy_Clinical/phs000915.rds")

# ====== phs001141 =======
phs001141 <- res_list$phs001141

phs001141 %<>%
  filter(IS_TUMOR != "normal", !is.na(SUBJECT_ID))
saveRDS(phs001141, file = "data/Tidy_Clinical/phs001141.rds")
