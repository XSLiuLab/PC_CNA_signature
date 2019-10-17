library(tidyverse)
library(readxl)
library(magrittr)

load("dbGap/dbGap_phenotype.RData")


# Clean datasets one by one -----------------------------------------------

# phs000447
phs000447 = res_list$phs000447
phs000447_T1 = read_excel(path = "data/TIdy_Clinical_Excel/phs000447.xlsx", sheet = 1)
phs000447_T2 = read_excel(path = "data/TIdy_Clinical_Excel/phs000447.xlsx", sheet = 2)
phs000447_T3 = read_excel(path = "data/TIdy_Clinical_Excel/phs000447.xlsx", sheet = 3)

phs000447_T1 = phs000447_T1 %$%
  str_remove_all(PatientID, "-Tumor") %>%
  str_remove_all("-Normal") %>%
  unique()

phs000447_T3 = phs000447_T3 %>%
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
