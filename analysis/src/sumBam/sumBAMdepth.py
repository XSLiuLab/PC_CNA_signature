#!/usr/bin/env python

# Author: Shixang Wang
# Copyright @ MIT
# Update date: 2020/03/07
import glob
from multiprocessing import Pool

DEPTH_DIRS = [
    "/public/home/liuxs/ncbi/dbGaP-16533/dnaseq/BQSR/bqsrbam/depth/",
    "/public/home/liuxs/ncbi/dbGaP-21926/dnaseq/BQSR/bqsrbam/depth/",
    "/public/home/liuxs/biodata/gdc/links/TCGA_PRAD/depth/"
]


def process_file(fp, reads_cutoff = 25):
    count = 0
    s = 0
    print("==> Processing file %s" % fp)
    with open(fp, "r", encoding="utf-8") as f:
        for row in f:
            content = row.split("\t")
            reads = int(content[2])
            if reads > reads_cutoff:
                count += 1
                s += reads
    print("file:", fp)
    print("total counts:", count)
    print("total reads:", s)
    if count < 1e6:
        if reads_cutoff <= 5:
            print("==> Still got few counts, return it.")
            write_result(fp, count, s)
            return True
        else:
            print("==> Too few counts, assume encount a normal sample, re-counting with lower cutoff...")
            return process_file(fp, reads_cutoff = 2)
    else:
        write_result(fp, count, s)
        return True

def write_result(fp, count, s):
    with open("/public/home/liuxs/PRAD_PROJ_BAM_SUMMARY.csv", "a") as f:
        f.write(fp + "," + str(count) + "," + str(s) + "\n")

def main():

    ## Create list for storing the results
    path_list = []

    ## Get all files
    for DIR in DEPTH_DIRS:
        print("=> Handling directory %s" % DIR)
        path_list += glob.glob(DIR + "/*depth")

    print("Total", len(path_list), "files found.")

    ## Use map
    with Pool(20) as p:
        p.map(process_file, path_list)

if __name__ == "__main__":
    main()