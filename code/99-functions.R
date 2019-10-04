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


