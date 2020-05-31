#! /usr/bin/env Rscript

library("sequenza")

args <- commandArgs(TRUE)
sample_id <- args[1]
input_file <- args[2]
out_dir <- args[3]

# https://github.com/ShixiangWang/copynumber
test <- sequenza.extract(input_file, assembly = "hg38")

CP <- sequenza.fit(test, female = FALSE)

sequenza.results(
  sequenza.extract = test,
  cp.table = CP,
  sample.id = sample_id,
  out.dir = out_dir,
  female = FALSE
)
