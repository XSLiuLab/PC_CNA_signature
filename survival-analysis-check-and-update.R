library(ezcox)
library(ggplot2)

df.seqz = readRDS(file = "output/df.seqz.RDS")
surv_df <- readRDS(file = "data/PRAD_Survival.rds")

range01 <- function(x, ...) {
  (x - min(x, ...)) / (max(x, ...) - min(x, ...))
}
cols_to_sigs.seqz <- c(paste0("CN-Sig", 1:5), paste0("SBS-Sig", 1:3))
surv_dt <- dplyr::inner_join(
  df.seqz %>%
    dplyr::select(-c("Study", "subject_id",
                     "tumor_body_site",
                     "tumor_Run", "normal_Run",
                     "CNV_ID", "Fusion")),
  surv_df, by = c("PatientID" = "sample"))
dup_ids = which(duplicated(surv_dt$PatientID))
# Remove duplicated records
surv_dt = surv_dt[-dup_ids, ]
# Scale the signature exposure to 0-20.
# To better understand the effect of signature exposure on clinical events,
# here we scale them into range 0-20. 1 increase means 5% increase of signature exposure.
#
# Scale the CNA burden to 0-20.
#
# keep only ti_fraction, tv result will be 1/ti_result
# also scale it to 0-20
# also scale purity to 0-20
surv_dt = surv_dt %>%
  dplyr::select(-Tv_fraction) %>%
  dplyr::mutate(cnaBurden = 20 * cnaBurden,
                Stage = as.character(Stage) %>%
                  factor(levels = c("T2", "T3", "T4"))) %>%
  dplyr::mutate_at(c(cols_to_sigs.seqz, "Ti_fraction", "purity"),
                   ~ 20 * range01(., na.rm = TRUE))



# Unvariable analysis -----------------------------------------------------

p = show_forest(surv_dt, covariates = cols_to_sigs.seqz,
                time = "OS.time", status = "OS",
                merge_models = TRUE, add_caption = FALSE)
p

# split into primary and meta
table(surv_dt$sample_type)

p1 = show_forest(surv_dt %>%
                   subset(sample_type == "Primary"),
                 covariates = cols_to_sigs.seqz,
                time = "OS.time", status = "OS",
                merge_models = TRUE, add_caption = FALSE)
p1

ggsave(filename = "check/OS_sig_primary.pdf", plot = p1,
       width = 7, height = 5)

p2 = show_forest(surv_dt %>%
                   subset(sample_type != "Primary"),
                 covariates = cols_to_sigs.seqz,
                 time = "OS.time", status = "OS",
                 merge_models = TRUE, add_caption = FALSE)
p2
ggsave(filename = "check/OS_sig_meta.pdf", plot = p2,
       width = 7, height = 5)

# Multivariable analysis
p = show_forest(surv_dt,
                covariates = cols_to_sigs.seqz[1],
                controls = cols_to_sigs.seqz[-1],
                time = "OS.time", status = "OS",
                add_caption = FALSE, merge_models = TRUE)
p
ggsave(filename = "check/OS_sig_multivariable.pdf", plot = p,
       width = 7, height = 5)

p1 = show_forest(surv_dt %>%
                   subset(sample_type == "Primary"),
                covariates = cols_to_sigs.seqz[1],
                controls = cols_to_sigs.seqz[-1],
                time = "OS.time", status = "OS",
                add_caption = FALSE, merge_models = TRUE)
p1
ggsave(filename = "check/OS_sig_multivariable_primary.pdf", plot = p1,
       width = 7, height = 5)

p2 = show_forest(surv_dt %>%
                   subset(sample_type != "Primary"),
                 covariates = cols_to_sigs.seqz[1],
                 controls = cols_to_sigs.seqz[-1],
                 time = "OS.time", status = "OS",
                 add_caption = FALSE, merge_models = TRUE)
p2
ggsave(filename = "check/OS_sig_multivariable_meta.pdf", plot = p2,
       width = 7, height = 5)


# PFS
p = show_forest(surv_dt, covariates = cols_to_sigs.seqz,
                time = "PFI.time", status = "PFI",
                merge_models = TRUE, add_caption = FALSE)
p

p1 = show_forest(surv_dt,
                 covariates = cols_to_sigs.seqz[1],
                 controls = cols_to_sigs.seqz[-1],
                time = "PFI.time", status = "PFI",
                merge_models = TRUE, add_caption = FALSE)
p1

ggsave(filename = "check/PFS_sig_multivariable.pdf", plot = p1,
       width = 7, height = 5)



# Features ----------------------------------------------------------------

cols_to_features <- c(
  "Age", "Stage", "GleasonScore",
  "n_SBS",
  "n_INDEL",
  "n_CNV", "n_Amp", "n_Del",
  "cnaBurden",
  "Ti_fraction",
  "TDP score",
  "Chromoth_state",
  "MATH",
  "purity",
  "ploidy"
)
p = show_forest(surv_dt,
                covariates = cols_to_features,
                controls = NULL,
                time = "OS.time", status = "OS",
                merge_models = TRUE, limits = log(c(0.01, 5)),
                add_caption = FALSE)
p

ggsave(filename = "check/OS_features.pdf", plot = p,
       width = 8, height = 6)

p1 = show_forest(surv_dt %>%
                  subset(sample_type == "Primary"),
                covariates = cols_to_features,
                controls = NULL,
                time = "OS.time", status = "OS",
                merge_models = TRUE, limits = log(c(0.01, 5)),
                add_caption = FALSE)
p1

ggsave(filename = "check/OS_features_primary.pdf", plot = p1,
       width = 8, height = 6)

p2 = show_forest(surv_dt %>%
                   subset(sample_type != "Primary"),
                 covariates = cols_to_features,
                 controls = NULL,
                 time = "OS.time", status = "OS",
                 merge_models = TRUE, limits = log(c(0.01, 5)),
                 add_caption = FALSE)
p2

ggsave(filename = "check/OS_features_meta.pdf", plot = p2,
       width = 8, height = 6)

# multivariable
p3 = show_forest(surv_dt,
                 covariates = cols_to_features[1],
                 controls = cols_to_features[-1],
                 time = "OS.time", status = "OS",
                 merge_models = TRUE, limits = log(c(0.01, 5)),
                 add_caption = FALSE)
p3

ggsave(filename = "check/OS_features_multivariable.pdf", plot = p3,
       width = 8, height = 6)

# PFI
p = show_forest(surv_dt,
                covariates = cols_to_features,
                controls = NULL,
                time = "PFI.time", status = "PFI",
                merge_models = TRUE, limits = log(c(0.5, 5)),
                add_caption = FALSE)
p

ggsave(filename = "check/PFS_features.pdf", plot = p,
       width = 8, height = 6)
p1 = show_forest(surv_dt %>%
                   subset(sample_type == "Primary"),
                 covariates = cols_to_features,
                 controls = NULL,
                 time = "PFI.time", status = "PFI",
                 merge_models = TRUE, limits = log(c(0.5, 5)),
                 add_caption = FALSE)
p1

ggsave(filename = "check/PFI_features_primary.pdf", plot = p1,
       width = 8, height = 6)

p2 = show_forest(surv_dt %>%
                   subset(sample_type != "Primary"),
                 covariates = cols_to_features,
                 controls = NULL,
                 time = "OS.time", status = "OS",
                 merge_models = TRUE, limits = log(c(0.01, 5)),
                 add_caption = FALSE)
p2

ggsave(filename = "check/OS_features_meta.pdf", plot = p2,
       width = 8, height = 6)

# multivariable
p3 = show_forest(surv_dt,
                 covariates = cols_to_features[1],
                 controls = cols_to_features[-1],
                 time = "PFI.time", status = "PFI",
                 merge_models = TRUE, limits = log(c(0.5, 5)),
                 add_caption = FALSE)
p3

ggsave(filename = "check/PFI_features_multivariable.pdf", plot = p3,
       width = 8, height = 6)

