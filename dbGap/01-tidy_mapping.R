library(tidyverse)


# Read run table from dbGap and select columns -----

cols <- c(
  "AvgSpotLen", "InsertSize", "Platform", "Instrument", "LibraryLayout", "Run", "SRA_Sample",
  "SRA_Study", "Sample_Name", "body_site", "histological_type",
  "gap_accession", "submitted_subject_id", "is_tumor", "study_design", "study_disease"
)

run1 <- data.table::fread("dbGap/SraRunTable_1.txt")[Assay_Type == "WXS"] %>%
  as_tibble() %>%
  select_at(cols)
run2 <- data.table::fread("dbGap/SraRunTable_2.txt")[Assay_Type == "WXS"] %>%
  as_tibble() %>%
  select_at(cols) %>%
  mutate(submitted_subject_id = as.character(submitted_subject_id))

runInfo <- bind_rows(run1, run2)
rm(run1, run2)

# Total unique submitted patient id
runInfo %>%
  mutate(ID = paste0(gap_accession, submitted_subject_id)) %>%
  pull(ID) %>%
  unique() %>%
  length()

# Clean ID mappings from sample pairs to bam files -----------------------------

paired_df <- runInfo %>%
  group_by(gap_accession, submitted_subject_id) %>%
  summarise(hasPaired = all(c("Yes", "No") %in% is_tumor), count = n()) %>%
  arrange(desc(count))

paired_ids <- paired_df %>%
  filter(hasPaired) %>%
  select(gap_accession, submitted_subject_id)

gen_mapping <- function(data) {
  # 只保留一个正常样本作为 control
  # 肿瘤样本按 Body site （NA 应该就是指前列腺）去重后与 control 配对

  filter_normal <- function(data) {
    message("==> Filtering normal...")
    data <- filter(data, is_tumor == "No")
    data_keep <- data %>%
      filter(grepl("blood", Sample_Name, ignore.case = TRUE) | grepl("blood", body_site, ignore.case = TRUE))
    if (nrow(data_keep) != 0) {
      data <- data_keep
    }

    if (nrow(data) > 1) {
      data <- head(data, 1)
    }

    data
  }


  filter_sample <- function(data) {
    # data_bk = data
    data_normal <- filter(data, is_tumor == "No")
    if (nrow(data_normal) > 1) {
      data_normal <- filter_normal(data_normal)
    }
    data_tumor <- filter(data, is_tumor == "Yes")
    if (nrow(data_tumor) > 1) {
      message("==> Filtering tumor...")

      data_tumor <- data_tumor %>%
        arrange(SRA_Sample) %>% # order by SRA Sample ID firstly
        distinct(body_site, .keep_all = TRUE)
    }

    bind_rows(data_normal, data_tumor)
  }

  message("=> Processing ", unique(data$submitted_subject_id))


  if (nrow(data) > 2) {
    data <- data %>%
      filter_sample()
  }

  # Create mapping

  data <- data %>%
    filter(is_tumor == "Yes") %>%
    transmute(
      gap_accession = gap_accession,
      subject_id = submitted_subject_id,
      tumor_sample = Sample_Name,
      normal_sample = filter(data, is_tumor == "No") %>% pull(Sample_Name),
      tumor_body_site = body_site,
      normal_body_site = filter(data, is_tumor == "No") %>% pull(body_site),
      tumor_SRA_ID = SRA_Sample,
      normal_SRA_ID = filter(data, is_tumor == "No") %>% pull(SRA_Sample),
      tumor_Run = Run,
      normal_Run = filter(data, is_tumor == "No") %>% pull(Run)
    )

  data
}

mapping_df <- paired_ids %>%
  rowwise() %>%
  do(df = filter(
    runInfo, gap_accession == .$gap_accession,
    submitted_subject_id == .$submitted_subject_id,
    !is.na(is_tumor)
  ) %>%
    gen_mapping()) %>%
  pull(df) %>%
  bind_rows() %>%
  arrange(gap_accession, subject_id)

mapping_df %>%
  mutate(ID = paste0(gap_accession, subject_id)) %>%
  pull(ID) %>%
  unique() %>%
  length()

save(mapping_df, file = "dbGap/mapping_df.RData")
write_tsv(mapping_df, path = "dbGap/mapping_df.tsv")

# ** Find a special case '511503*', test it ---------
# 到底该保留 InsertSize 0 还是 ~220 左右的？

zz <- filter(runInfo, startsWith(submitted_subject_id, "511503")) %>%
  select(Run, SRA_Sample, submitted_subject_id, Sample_Name, is_tumor, body_site, Instrument, InsertSize)
# Source filter_* functions
zz

tibble(id = unique(zz$submitted_subject_id)) %>%
  rowwise() %>%
  do(df = filter(zz, submitted_subject_id == .$id) %>%
    filter_normal()) %>%
  pull(df) %>%
  bind_rows()

tibble(id = unique(zz$submitted_subject_id)) %>%
  rowwise() %>%
  do(df = filter(zz, submitted_subject_id == .$id) %>%
    filter_sample()) %>%
  pull(df) %>%
  bind_rows()

# Read phenotype data for subjects and samples ----------------------------
read_dbGap <- function(accession, destdir = getwd(),
                       type = c("subject", "sample"),
                       col_types = cols(.default = "c"),
                       pattern_subject_file = "Subject_Phenotypes",
                       pattern_sample_file = "Sample_Attributes") {
  # col_types see ?readr::read_tsv
  # user can set it to dplyr::cols() if no error 'Too many conversion specifiers in format string' returned

  type <- match.arg(type)
  # append / at the end if user not type it
  destdir <- ifelse(endsWith(destdir, suffix = "/"), destdir, paste0(destdir, "/"))

  if (type == "subject") {
    fl <- dir(destdir, paste0(accession, ".*", pattern_subject_file), full.names = TRUE)
  } else {
    fl <- dir(destdir, paste0(accession, ".*", pattern_sample_file), full.names = TRUE)
  }

  if (length(fl) < 1) stop("Cannot find any file! Please check your input.")

  message("==> Finding out files with type ", type)
  print(fl)

  if (!require(tidyverse)) {
    message("Please install tidyverse package and re-run this function!")
  } else {
    purrr::map(fl, ~ readr::read_tsv(.,
      comment = "#", progress = TRUE,
      col_types = col_types
    )) %>%
      dplyr::bind_rows() %>%
      unique()
  }
}

read_dbGapList <- function(accession_list, destdir = getwd(),
                           col_types = cols(.default = "c"),
                           pattern_subject_file = "Subject_Phenotypes",
                           pattern_sample_file = "Sample_Attributes") {
  if (!require(tidyverse)) {
    message("Please install tidyverse package and re-run this function!")
  } else {
    purrr::map(accession_list, function(x) {
      message("=> Processing accession ", x)
      type <- c("subject", "sample")
      purrr::map(
        type,
        ~ suppressWarnings(read_dbGap(x, destdir = destdir, type = .))
      ) %>%
        setNames(type)
    }) %>%
      setNames(accession_list)
  }
}

df_list <- read_dbGapList(
  accession_list = paste0(
    "phs00",
    c("0447", "0554", "0909", "0915", "1141")
  ),
  destdir = "dbGap/phenotype/"
)




# Clean phenotype data for studies one by one -----
ns <- names(df_list)

# ** phs000447 ---------
ns[1]
res_list <- list()

df_list$phs000447$subject %>% View()

tmp <- df_list$phs000447$sample
tmp <- tmp %>%
  select(-c("Original Material Type", "Sample Zonal Origin")) %>%
  mutate(SUBJID = sub("([^_]+)_.*", "\\1", SAMPID)) %>% # Remove _*
  mutate(SUBJID = ifelse(startsWith(SUBJID, "STID"),
    SUBJID,
    gsub("[^0-9-]", "", SUBJID)
  )) # Remove any non-numeric characters

res_list[[ns[1]]] <- df_list$phs000447$subject %>%
  full_join(tmp, by = "SUBJID")

res_list$phs000447 %>% View()

# ** phs000554 ---------
ns[2]

df_list$phs000554$subject %>% View()

tmp <- df_list$phs000554$sample

bridge <- runInfo %>%
  filter(gap_accession == ns[2]) %>%
  select(submitted_subject_id, Sample_Name) %>%
  unique()

res_list[[ns[2]]] <- df_list$phs000554$subject %>%
  full_join(bridge, by = c("SUBJECT_ID" = "submitted_subject_id")) %>%
  full_join(tmp, by = c("Sample_Name" = "SAMPLE_ID"))


res_list$phs000554 %>% View()


# ** phs000909 ---------
ns[3]

df_list$phs000909$subject %>% View()

tmp <- df_list$phs000909$sample
tmp <- tmp %>%
  select(-c("IN_NAT_MED_PAPER")) %>%
  mutate(SUBJECT_ID = sub("([^_]+)_.*", "\\1", SAMPLE_ID))

res_list[[ns[3]]] <- df_list$phs000909$subject %>%
  full_join(tmp, by = "SUBJECT_ID")

res_list$phs000909 %>% View()

# ** phs000915 ---------
ns[4]

df_list$phs000915$subject %>% View()

tmp <- df_list$phs000915$sample

bridge <- runInfo %>%
  filter(gap_accession == ns[4]) %>%
  select(submitted_subject_id, Sample_Name) %>%
  unique()

"1115154" %in% bridge$submitted_subject_id


res_list[[ns[4]]] <- df_list$phs000915$subject %>%
  mutate(SUBJECT_ID = sub("^0", "", SUBJECT_ID)) %>%
  full_join(bridge, by = c("SUBJECT_ID" = "submitted_subject_id")) %>%
  full_join(tmp, by = c("Sample_Name" = "SAMPLE_ID"))

res_list$phs000915 %>% View()

# ** phs001141 ---------
ns[5]

df_list$phs001141$subject %>% View()

tmp <- df_list$phs001141$sample

bridge <- runInfo %>%
  filter(gap_accession == ns[5]) %>%
  select(submitted_subject_id, Sample_Name) %>%
  unique()


res_list[[ns[5]]] <- df_list$phs001141$subject %>%
  full_join(bridge, by = c("SUBJECT_ID" = "submitted_subject_id")) %>%
  full_join(tmp, by = c("Sample_Name" = "SAMPLE_ID"))

res_list$phs001141 %>% View()


# Save data ---------------------------------------------------------------

save(res_list, file = "dbGap/dbGap_phenotype.RData")
