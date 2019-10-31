#// Homebrew functions to (pre)process data

# FUN: extract copy number values from FACETS results ---------------------

#// The commands in this function are mostly created by Minfang
#// and then packaged by Shixiang
# extract_facets_cnv = function(target_dir, target_path) {
#   SAMPLE = dir(target_dir, pattern = "fit.RData") %>%
#     str_remove(".fit.RData")
#
#   cols <- c("chrom","start","end","tcn.em","lcn.em")
#   all_sample <- rbind()
#   for (sample in SAMPLE){
#     sample_RData <- file.path(target_dir, paste(sample, ".fit.RData", sep = ""))
#     load(sample_RData)
#     df <- fit$cncf[cols]
#     df["sample"] <- sample
#     all_sample <- rbind(all_sample, df)
#     warnings = paste("processing ",sample,sep = "")
#     message(warnings)
#   }
#
#   colnames(all_sample) <- c("Chromosome","Start.bp","End.bp","modal_cn","minor_cn","sample")
#   facets_CNV <- all_sample %>% select(Chromosome, Start.bp,End.bp,modal_cn,sample)
#   write.table(facets_CNV, target_path, sep = "\t", quote = FALSE, row.names = F)
#   rm(list = ls())
# }

# Process data from HuiMin
extract_facets_cnv = function(target_dir, target_path) {
  SAMPLE = dir(target_dir, pattern = ".Rdata") %>%
    str_remove(".Rdata")

  cols <- c("chrom","start","end","tcn.em","lcn.em")
  all_sample <- rbind()
  for (sample in SAMPLE){
    sample_RData <- file.path(target_dir, paste(sample, ".Rdata", sep = ""))
    load(sample_RData)
    df <- fit$cncf[cols]
    df["sample"] <- sample
    all_sample <- rbind(all_sample, df)
    warnings = paste("processing ",sample,sep = "")
    message(warnings)
  }

  colnames(all_sample) <- c("Chromosome","Start.bp","End.bp","modal_cn","minor_cn","sample")
  facets_CNV <- all_sample %>% select(Chromosome, Start.bp,End.bp,modal_cn,sample)
  write.table(facets_CNV, target_path, sep = "\t", quote = FALSE, row.names = F)
  rm(list = ls())
}

extract_facets_purity_and_ploidy = function(target_dir, target_path) {
  SAMPLE = dir(target_dir, pattern = ".Rdata") %>%
    str_remove(".Rdata")

  all_sample <- rbind()
  for (sample in SAMPLE){
    sample_RData <- file.path(target_dir, paste(sample, ".Rdata", sep = ""))
    load(sample_RData)
    df = data.frame(
      sample = sample,
      purity = round(as.numeric(fit$purity), 2),
      ploidy = round(as.numeric(fit$ploidy), 2),
      stringsAsFactors = FALSE)
    all_sample <- rbind(all_sample, df)
    warnings = paste("processing ",sample,sep = "")
    message(warnings)
  }

  colnames(all_sample) <- c("sample", "purity", "ploidy")
  write.table(all_sample, target_path, sep = "\t", quote = FALSE, row.names = F)
  rm(list = ls())
}

facets_to_GISTIC2 = function(target_dir, target_path) {
  SAMPLE = dir(target_dir, pattern = ".Rdata") %>%
    str_remove(".Rdata")

  cols <- c("chrom","start","end","num.mark","tcn.em")
  all_sample <- rbind()
  for (sample in SAMPLE){
    sample_RData <- file.path(target_dir, paste(sample, ".Rdata", sep = ""))
    load(sample_RData)
    df <- fit$cncf[cols]
    df["sample"] <- sample
    all_sample <- rbind(all_sample, df)
    warnings = paste("processing ",sample,sep = "")
    message(warnings)
  }

  colnames(all_sample) <- c("chrom","loc.start","loc.end", "num.mark", "seg.mean","ID")
  facets_CNV <- all_sample %>% dplyr::select(ID, dplyr::everything())
  facets_CNV$seg.mean = ifelse(facets_CNV$seg.mean!=0,
                               log2(facets_CNV$seg.mean) - 1, -11) # set a default minimum seg.mean for which copy number equals 0
  # -11 comes from PRAD 1000 study minumum seg.mean
  write.table(facets_CNV, target_path, sep = "\t", quote = FALSE, row.names = F)
  rm(list = ls())
}

# FUN: visualize drake plan quickly ---------------------------------------

vis_plan = function(plan,
                    file = character(0),
                    selfcontained = FALSE,
                    digits = 1,
                    targets_only = FALSE,
                    font_size = 15,
                    main = NULL, ...) {
  config = drake_config(plan)
  vis_drake_graph(config,
                  file = file,
                  selfcontained = selfcontained,
                  digits = digits,
                  targets_only = targets_only,
                  font_size = font_size,
                  main = main, ...)
}


