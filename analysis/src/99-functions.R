# // Homebrew functions to (pre)process data

# FUN: extract copy number values from FACETS results ---------------------

# // The commands in this function are mostly created by Minfang
# // and then packaged by Shixiang
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
extract_facets_cnv <- function(target_dir, target_path) {
  SAMPLE <- dir(target_dir, pattern = ".Rdata") %>%
    str_remove(".Rdata")

  cols <- c("chrom", "start", "end", "tcn.em", "lcn.em")
  all_sample <- rbind()
  for (sample in SAMPLE) {
    sample_RData <- file.path(target_dir, paste(sample, ".Rdata", sep = ""))
    load(sample_RData)
    df <- fit$cncf[cols]
    df["sample"] <- sample
    all_sample <- rbind(all_sample, df)
    warnings <- paste("processing ", sample, sep = "")
    message(warnings)
  }

  colnames(all_sample) <- c("Chromosome", "Start.bp", "End.bp", "modal_cn", "minor_cn", "sample")
  facets_CNV <- all_sample %>% dplyr::select(Chromosome, Start.bp, End.bp, modal_cn, sample)
  write.table(facets_CNV, target_path, sep = "\t", quote = FALSE, row.names = F)
  rm(list = ls())
}

# FUN: extract purity and ploidy from FACETS results ---------------------

extract_facets_purity_and_ploidy <- function(target_dir, target_path) {
  SAMPLE <- dir(target_dir, pattern = ".Rdata") %>%
    str_remove(".Rdata")

  all_sample <- rbind()
  for (sample in SAMPLE) {
    sample_RData <- file.path(target_dir, paste(sample, ".Rdata", sep = ""))
    load(sample_RData)
    df <- data.frame(
      sample = sample,
      purity = round(as.numeric(fit$purity), 2),
      ploidy = round(as.numeric(fit$ploidy), 2),
      stringsAsFactors = FALSE
    )
    all_sample <- rbind(all_sample, df)
    warnings <- paste("processing ", sample, sep = "")
    message(warnings)
  }

  colnames(all_sample) <- c("sample", "purity", "ploidy")
  write.table(all_sample, target_path, sep = "\t", quote = FALSE, row.names = F)
  rm(list = ls())
}

# FUN: generate GISTIC2 input data from FACETS results ---------------------
# See https://github.com/ShixiangWang/install_GISTIC to install GISTIC2

facets_to_GISTIC2 <- function(target_dir, target_path, rm_X = TRUE, rm_samps = NULL) {
  SAMPLE <- dir(target_dir, pattern = ".Rdata") %>%
    str_remove(".Rdata")

  cols <- c("chrom", "start", "end", "num.mark", "tcn.em")
  all_sample <- rbind()
  for (sample in SAMPLE) {
    sample_RData <- file.path(target_dir, paste(sample, ".Rdata", sep = ""))
    load(sample_RData)
    df <- fit$cncf[cols]
    df["sample"] <- sample
    all_sample <- rbind(all_sample, df)
    warnings <- paste("processing ", sample, sep = "")
    message(warnings)
  }

  colnames(all_sample) <- c("chrom", "loc.start", "loc.end", "num.mark", "seg.mean", "ID")
  facets_CNV <- all_sample %>% dplyr::select(ID, dplyr::everything())
  facets_CNV$seg.mean <- ifelse(facets_CNV$seg.mean != 0,
    log2(facets_CNV$seg.mean) - 1, -11
  ) # set a default minimum seg.mean for which copy number equals 0
  # -11 comes from PRAD 1000 study minumum seg.mean

  # Remove chr23 (chrX) due to patients are male
  # 1 copy will be treated as deletions
  if (rm_X) {
    facets_CNV <- subset(facets_CNV, chrom != 23)
  }

  if (!is.null(rm_samps)) {
    facets_CNV <- subset(facets_CNV, !ID %in% rm_samps)
  }

  write.table(facets_CNV, target_path, sep = "\t", quote = FALSE, row.names = F)
  rm(list = ls())
}


# FUN: extract copy number values from Sequenza results ---------------------

extract_seqz_cnv <- function(target_dir, target_path) {
  SAMPLE <- dir(target_dir, pattern = "_segments.txt") %>%
    stringr::str_remove("_segments.txt")
  res <- purrr::map2_df(file.path(
    target_dir,
    paste(SAMPLE,
          "_segments.txt",
          sep = ""
    )
  ), SAMPLE, function(x, y) {
    message("Processing ", y)
    df <- readr::read_tsv(x)
    df <- df %>%
      dplyr::select(chromosome, start.pos, end.pos, CNt) %>%
      dplyr::mutate(sample = y)
    colnames(df) <- c("Chromosome", "Start.bp", "End.bp", "modal_cn", "sample")
    df
  })
  readr::write_tsv(res, path = target_path)
}


# FUN: transform sequenza data to GISTIC2 input format --------------------


seqz_to_GISTIC2 <- function(target_dir, target_path, rm_X = TRUE, rm_samps = NULL) {
  SAMPLE <- dir(target_dir, pattern = "_segments.txt") %>%
    stringr::str_remove("_segments.txt")
  res <- purrr::map2_df(file.path(
    target_dir,
    paste(SAMPLE,
          "_segments.txt",
          sep = ""
    )
  ), SAMPLE, function(x, y) {
    message("Processing ", y)
    df <- readr::read_tsv(x, col_types = readr::cols())
    df <- df %>%
      dplyr::select(chromosome, start.pos, end.pos, N.BAF, CNt) %>%
      dplyr::mutate(sample = y)
    colnames(df) <- c("chrom", "loc.start", "loc.end", "num.mark", "seg.mean", "ID")
    df
  })

  seqz_CNV <- res %>% dplyr::select(ID, dplyr::everything())
  seqz_CNV$seg.mean <- ifelse(seqz_CNV$seg.mean != 0,
                              log2(seqz_CNV$seg.mean) - 1, -11
  ) # set a default minimum seg.mean for which copy number equals 0
  # -11 comes from PRAD 1000 study minumum seg.mean

  # Remove chr23 (chrX) due to patients are male
  # 1 copy will be treated as deletions
  if (rm_X) {
    seqz_CNV <- subset(seqz_CNV, chrom != 23)
  }

  if (!is.null(rm_samps)) {
    seqz_CNV <- subset(seqz_CNV, !ID %in% rm_samps)
  }

  write.table(seqz_CNV, target_path, sep = "\t", quote = FALSE, row.names = F)
}


# FUN: extract purity and ploidy from Sequenza results ---------------------

extract_seqz_purity_and_ploidy <- function(target_dir, target_path) {
  SAMPLE <- dir(target_dir, pattern = "_alternative_solutions.txt") %>%
    stringr::str_remove("_alternative_solutions.txt")
  res <- purrr::map2_df(file.path(
    target_dir,
    paste(SAMPLE,
          "_alternative_solutions.txt",
          sep = ""
    )
  ), SAMPLE, function(x, y) {
    message("Processing ", y)
    df <- readr::read_tsv(x)
    df <- df %>%
      dplyr::select(cellularity, ploidy) %>%
      dplyr::mutate(sample = y) %>%
      dplyr::slice(1) # Choose the best solution
    colnames(df) <- c("purity", "ploidy", "sample")
    df
  })
  readr::write_tsv(res, path = target_path)
}


# FUN: open directory ---------------------
# from rvcheck package
o <- function(file = ".") {
  is_r_server <- function() {
    if (!exists("RStudio.Version")) {
      return(FALSE)
    }
    if (!is.function(RStudio.Version)) {
      return(FALSE)
    }
    RStudio.Version()$mode == "server"
  }

  if (is_r_server()) {
    while (Sys.getenv("rserver_ip") == "") {
      Sys.setenv(rserver_ip = readline("Set your RStudio Server IP address: "))
    }
    cat("Got!")
    while (Sys.getenv("rserver_port") == "") {
      Sys.setenv(rserver_port = readline("Set your RStudio Server port: "))
    }
    cat("Got!")

    if (!file.exists(file)) {
      stop("File ", file, " does not exist!")
    } else {
      rserver_ip <- Sys.getenv("rserver_ip")
      rserver_port <- Sys.getenv("rserver_port")
      if (!startsWith(rserver_ip, "http")) {
        rserver_ip <- paste0("http://", rserver_ip)
      }
      utils::browseURL(
        paste0(
          paste(rserver_ip, rserver_port, sep = ":"),
          "/file_show?path=",
          file
        )
      )
    }
  } else {
    os <- Sys.info()[1]
    if (os == "Darwin") {
      cmd <- paste("open", file)
      system(cmd)
    } else if (os == "Linux") {
      cmd <- paste("xdg-open", file, "&")
      system(cmd)
    } else if (os == "Windows") {
      ## wd <- sub("/", "\\", getwd())
      ## cmd <- paste("explorer", wd)
      ## suppressWarnings(shell(cmd))
      cmd <- paste("start", file)
      shell(cmd)
    }
  }
}
