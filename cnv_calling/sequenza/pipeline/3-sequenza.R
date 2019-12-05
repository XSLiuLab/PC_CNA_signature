#! /usr/bin/env Rscript

library("sequenza")

args = commandArgs(TRUE)
sample_id = args[1]
input_file = args[2]
out_dir = args[3]

# https://github.com/ShixiangWang/copynumber
test = sequenza.extract(input_file, assembly = "hg38")

CP = sequenza.fit(test, female = FALSE)
#CP = sequenza.fit(test,segment.filter = 5e6, N.ratio.filter = 80, N.BAF.filter = 50) 

sequenza.results(
    sequenza.extract = test,
    cp.table = CP,
    sample.id = sample_id,
    out.dir = out_dir,
    female = FALSE)

