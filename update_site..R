#!/usr/bin/env Rscript

output_dir = "docs"

## Check the workding directory

if (!dir.exists("analysis")) {
  stop("Please check the working directory!")
} else {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  ## FILE PATH
  index_html = "analysis/index.html"
  index_files = "analysis/index_files"
  fig_files = c("analysis/fig/")#, "analysis/corrr_network.pdf")

  file.copy(c(index_html, index_files, fig_files), to = output_dir, recursive = TRUE)
  #file.remove(index_html)
  unlink(x = c(index_html, index_files, "analysis/index_cache"), recursive = TRUE)
}
