# Combine and clean IDs deeply
# using IDs from:
# 1) FACETS result or tidy mapping
# 2) 02-tidy_clinical_data.R
# 3) clinical data from Nat.Gen 2018 (1013 samples)
library(tidyverse)
library(sigminer)

# 1)
CNV = read_copynumber("data/CNV_from_dbGAP_PLUS_TCGA_WES_CVAL150.tsv", genome_build = "hg38",
                      complement = FALSE, verbose = TRUE)
FACETS_SAMPLES = CNV@summary.per.sample$sample
load("dbGap/mapping_df.RData")

# 2)
phs000447 = readRDS("data/Tidy_Clinical/phs000447.rds")
phs000554 = readRDS("data/Tidy_Clinical/phs000554.rds")
phs000909 = readRDS("data/Tidy_Clinical/phs000909.rds")
phs000915 = readRDS("data/Tidy_Clinical/phs000915.rds")
phs001141 = readRDS("data/Tidy_Clinical/phs001141.rds")

phs000554 = phs000554 %>%
  slice(1:61) %>%
  mutate(SampleName = sub(" (.*)", "", SampleName))

# 3)
NatGen_clincal <- readxl::read_excel("data/NIHMS938028-supplement-5.xlsx", skip = 2)
names(NatGen_clincal)[2] <- "Tumor_Sample_Barcode"

# NatGen_clincal = NatGen_clincal %>%
#   filter(Data.Source != "AAPC")

# Explore relationships ---------------------------------------------------
sum(phs000447$Sample %in% NatGen_clincal$Tumor_Sample_Barcode)
sum(phs000909$PatientID %in% NatGen_clincal$Tumor_Sample_Barcode) # strange
sum(phs000915$`PATIENT ID` %in% NatGen_clincal$Tumor_Sample_Barcode)

# It is very hard to merge all these datasets
# We focus the samples with raw sequence data processed by us


# Create unified clinical info --------------------------------------------

mapping_df %>%
  filter(tumor_body_site == "") %>%
  pull(gap_accession) %>%
  table()

mapping_df %>%
  filter(tumor_body_site == "N.A.") %>%
  pull(gap_accession) %>%
  table()

clean_maps = mapping_df %>%
  mutate(
    sample_type = case_when(
      tumor_body_site == "" ~ "Metastatic",
      tumor_body_site == "N.A." ~ "Unknown",  # Check it using other annotation data
      grepl("prostate", tumor_body_site, ignore.case = TRUE) ~ "Primary",
      TRUE ~ "Metastatic"
    )
  ) %>%
  select(gap_accession, subject_id, sample_type, tumor_body_site, tumor_Run, normal_Run)

table(clean_maps$sample_type)

# Add annotation one by one
#==== phs000447 ========
samps = clean_maps %>% filter(gap_accession == "phs000447") %>% pull(subject_id) %>% unique()
all(samps %in% phs000447$Predict_SampleID)

phs000447_T = phs000447 %>%
  select(-c(`dbGaP SubjID`, `Primary Disease`)) %>%
  rename(PSA = `PSA Pre-Operative`, RadiationTherapy=`Radiation Therapy`, Stage=`Cancer Stage`,
         GleasonScore = `Gleason Score`, Fusion = `TMPRSS2-ERG Fusion Status`) %>%
  mutate(Fusion = case_when(
    Fusion == "---" ~ NA_character_,
    grepl("NRF1-BRAF", Fusion) ~ "NRF1-BRAF",
    Fusion == "Negative" ~ "Negative",
    TRUE ~ "ERG"
  ))

# Found six replicate samples
# phs000447_T %>% select(Predict_SampleID) %>% duplicated() %>% which()

table(phs000447_T$Fusion)

#==== phs000554 ========
samps = clean_maps %>% filter(gap_accession == "phs000554") %>% pull(subject_id) %>% unique()
all(samps %in% phs000554$SUBJECT_ID)
samps[!samps %in% phs000554$SUBJECT_ID]

phs000554_T = phs000554 %>%
  mutate(SUBJECT_ID = ifelse(is.na(SUBJECT_ID), "44", SUBJECT_ID),
         PRIMARY_METASTATIC_TUMOR = ifelse(is.na(PRIMARY_METASTATIC_TUMOR), "Metastatic", PRIMARY_METASTATIC_TUMOR),
         SampleName = ifelse(SampleName == "WA43-44", "WA43", SampleName)) %>%
  select(-`Tumor/Normal location`) %>%
  filter(!SampleName %in% c("WA43-27", "WA43-71")) %>%
  rename(Fusion = `ETS/RAF/SPINK1 status`, PSA = `Serum PSA`,
         GleasonScore = `Gleason score`,
         PriorTreatment = `Prior Treatment`) %>%
  mutate(Fusion = case_when(
    grepl("ERG+", Fusion) ~ "ERG",
    grepl("ETV1+", Fusion) ~ "ERG",  # here ERG means ETS family, change it to ETV latter
    grepl("No ETS", Fusion) ~ "Negative",
    grepl("RAF1+", Fusion) ~ "RAF",
    grepl("SPINK1+", Fusion) ~ "SPINK1"
  ))

#==== phs000909 ========
samps = clean_maps %>% filter(gap_accession == "phs000909") %>% pull(subject_id) %>% unique()
all(samps %in% phs000909$PatientID)
samps[!samps %in% phs000909$PatientID]

phs000909_T = phs000909 %>%
  select(PatientID, `Pathology Classification`, `Purity (CLONET)`, `Ploidy (CLONET)`, Genomic_Burden, `Tumor RNASeq ID`, `Methylation ID`)

#==== phs000915 ========
samps = clean_maps %>% filter(gap_accession == "phs000915") %>% pull(subject_id) %>% unique()
all(samps %in% phs000915$`cBio_SU2C ID`)

phs000915_T = phs000915 %>%
  select(-`FIGURE 2 CASE#`, -`BIOPSY SITE`, -`SEQUENCING INSTITUTE`) %>%
  rename(SU2C_ID=`cBio_SU2C ID`, PatientID=`PATIENT ID`, Age=AGE,
         PRIOR_ABI_or_ENZ = `PRIOR ABI or ENZ`, PRIOR_TAXAN = `PRIOR TAXAN`,
         ClINICAL_SITE = `ClINICAL  SITE`, Purity = `TUMOR CONTENT ESTIMATED BY SEQ`)

#==== phs000915 ========
samps = clean_maps %>% filter(gap_accession == "phs001141") %>% pull(subject_id) %>% unique()
all(samps %in% phs001141$SUBJECT_ID)

phs001141_T = phs001141 %>%
  select(SUBJECT_ID, age, race, composite_progression, sex, Sample_Name, PRIMARY_METASTATIC_TUMOR,
         Stage, TUMOR_TREATMENT) %>%
  rename(Age = age, Race = race, Gender = sex)

save(phs000447_T, phs000554_T, phs000909_T, phs000915_T, phs001141_T, file = "data/dbGap_clean_phenotype.RData")


# Merge data --------------------------------------------------------------

# Merge all phenotype data from dbGap

dbGap_data =
  list(
    phs000447 = filter(clean_maps, gap_accession=="phs000447") %>%
      full_join(phs000447_T, by = c("subject_id"="Predict_SampleID")) %>%
      select(-Gender, -RadiationTherapy) %>%
      rename(PatientID = Sample),
    phs000554 = filter(clean_maps, gap_accession=="phs000554") %>%
      full_join(phs000554_T, by = c("subject_id"="SUBJECT_ID")) %>%
      rename(PatientID = SampleName) %>%
      select(gap_accession:PSA, PRIMARY_METASTATIC_TUMOR, -DiseaseState) %>%
      mutate(sample_type = PRIMARY_METASTATIC_TUMOR) %>% # correct sample_type
      select(-PRIMARY_METASTATIC_TUMOR),
    phs000909 = filter(clean_maps, gap_accession=="phs000909") %>%
      full_join(phs000909_T, by = c("subject_id"="PatientID")) %>%
      select(gap_accession:normal_Run),
    phs000915 = filter(clean_maps, gap_accession=="phs000915") %>%
      full_join(phs000915_T, by = c("subject_id"="SU2C_ID")) %>%
      select(gap_accession:Age) %>%
      mutate(Age = as.integer(Age)),
    phs001141 = filter(clean_maps, gap_accession=="phs001141") %>%
      full_join(phs001141_T, by = c("subject_id"="SUBJECT_ID")) %>%
      mutate(sample_type = PRIMARY_METASTATIC_TUMOR,
             Age = as.integer(Age)) %>%
      select(gap_accession:Age)
  ) %>% purrr::map(unique)
# Then merge data from Nat.Gen
# lapply(dbGap_data, colnames)
# sapply(dbGap_data, nrow)
# sapply(dbGap_data, function(x) x %>% unique() %>% nrow())
# sapply(dbGap_data, function(x) x %>% unique() %>% nrow()) %>% sum()
#purrr::reduce(dbGap_data, bind_rows)
dbGap = purrr::map_df(dbGap_data, bind_rows) %>%
  rename(Study = gap_accession)
table(dbGap$sample_type)
table(dbGap$Stage)
table(dbGap$Fusion)

# https://stackoverflow.com/questions/34096162/dplyr-mutate-replace-several-columns-on-a-subset-of-rows
mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}


# https://en.wikipedia.org/wiki/ETS_transcription_factor_family
dbGap %>%
  mutate_cond(grepl("RAF", Fusion), Fusion = "BRAF") %>%
  mutate_cond(grepl("ERG", Fusion), Fusion = "ETS") %>%
  filter(!is.na(Study)) %>%
  mutate(PatientID = ifelse(
    is.na(PatientID),
    paste(Study, subject_id, sep = "-"),
    PatientID
  ),
  GleasonScore = ifelse(grepl("N", GleasonScore), NA_integer_, GleasonScore),
  GleasonScore = ifelse(nchar(GleasonScore) >= 3, substr(GleasonScore, 1, 3), GleasonScore),
  GleasonScore = sapply(GleasonScore, function(x) eval(parse(text=x))) %>%
    as.integer()) -> zz

nrow(zz)
nrow(unique(zz))
table(clean_maps$gap_accession)
table(zz$Study)
