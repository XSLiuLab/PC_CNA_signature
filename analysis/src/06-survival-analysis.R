library(ezcox)

range01 <- function(x, ...) {
  (x - min(x, ...)) / (max(x, ...) - min(x, ...))
}

df.seqz = readRDS(file = "output/PRAD_Merge_Info_CNV_from_sequenza.RData")
surv_df <- readRDS(file = "data/PRAD_Survival.rds")

cols_to_sigs.seqz <- c(paste0("CNV_Sig", 1:5), paste0("SNV_Sig", 1:6))
table(surv_df$Study)

df.seqz$PatientID
surv_df$sample

surv_dt <- dplyr::inner_join(df.seqz %>%
                               dplyr::select(-c("Study", "subject_id", "tumor_body_site",
                                                "tumor_Run", "normal_Run", "CNV_ID",
                                                "Fusion", "sample_type")),
                             surv_df, by = c("PatientID" = "sample"))

dup_ids = which(duplicated(surv_dt$PatientID))
# Remove duplicated records
surv_dt = surv_dt[-dup_ids, ]

# Scale the signature activity to 0-20.
# To better understand the effect of signature activity on clinical events,
# here we scale them into range 0-20. 1 increase means 5% increase of signature activity.
#
# Scale the CNA burden to 0-20.
surv_dt = surv_dt %>%
  dplyr::mutate(cna_burden = 20 * cna_burden) %>%
  dplyr::mutate_at(cols_to_sigs.seqz, ~ 20 * range01(., na.rm = TRUE))

saveRDS(surv_dt, file = "output/PRAD_merged_survival_dt_seqz.rds")

# Unvariable analysis

show_forest(surv_dt, covariates = cols_to_sigs.seqz,
            time = "OS.time", status = "OS", merge_models = TRUE)

show_forest(surv_dt, covariates = cols_to_sigs.seqz,
            time = "PFI.time", status = "PFI", merge_models = TRUE)

show_forest(surv_dt %>% dplyr::filter(Study == "TCGA"), covariates = cols_to_sigs.seqz,
            time = "OS.time", status = "OS", merge_models = TRUE)

show_forest(surv_dt %>% dplyr::filter(Study == "phs000554"), covariates = cols_to_sigs.seqz,
            time = "OS.time", status = "OS", merge_models = TRUE)

# Multi-variable analysis

show_forest(surv_dt %>% dplyr::filter(Study == "TCGA"),
            covariates = cols_to_sigs.seqz, controls = c("purity", "Stage"),
            time = "OS.time", status = "OS", merge_models = TRUE, drop_controls = TRUE)

show_forest(surv_dt %>% dplyr::filter(Study == "TCGA"),
            covariates = cols_to_sigs.seqz, controls = c("purity", "Age"),
            time = "OS.time", status = "OS", merge_models = TRUE, drop_controls = TRUE)

show_forest(surv_dt %>% dplyr::filter(Study == "TCGA"),
            covariates = cols_to_sigs.seqz, controls = c("purity", "Age"),
            time = "PFI.time", status = "PFI", merge_models = TRUE, drop_controls = TRUE)

show_forest(surv_dt %>% dplyr::filter(Study == "TCGA"),
            covariates = cols_to_sigs.seqz, controls = c("purity", "Stage"),
            time = "PFI.time", status = "PFI", merge_models = TRUE, drop_controls = TRUE)
