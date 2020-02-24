
# Call absolute copy number from seg data ---------------------------------

library(readr)
library(data.table)
library(DoAbsolute)


data_dir <- "/public/data/cbioPortal/"

#// prad su2c
seg_tumor <- file.path(data_dir, "prad_su2c_2019/data_cna_hg19.seg")            # segment file
maf_tumor <- file.path(data_dir, "prad_su2c_2019/data_mutations_extended.txt")  # MAF file

Seg <- fread(seg_tumor)
colnames(Seg)
colnames(Seg) <- c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")
Maf <- fread(maf_tumor)

# Run DoAbsolute
DoAbsolute(
  Seg = Seg, Maf = Maf, platform = "Illumina_WES",
  copy.num.type = "total",
  primary.disease = "PRAD",
  max.as.seg.count = 10000, # Keep all samples processed by ABSOLUTE
  temp.dir = "output/DoAbsolute_prad_su2c",
  results.dir = "output/DoAbsolute_prad_su2c", nThread = 10,
  keepAllResult = TRUE, verbose = TRUE
)


#// mskcc 2020
seg_tumor <- file.path(data_dir, "prad_mskcc_2020/data_cna_hg19.seg")            # segment file
maf_tumor <- file.path(data_dir, "prad_mskcc_2020/data_mutations_extended.txt")  # MAF file

Seg <- fread(seg_tumor)
colnames(Seg)
colnames(Seg) <- c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")
Maf <- fread(maf_tumor, skip = 2)
colnames(Maf)

# Run DoAbsolute
DoAbsolute(
  Seg = Seg, Maf = Maf, platform = "Illumina_WES",
  copy.num.type = "total",
  primary.disease = "PRAD",
  max.as.seg.count = 10000, # Keep all samples processed by ABSOLUTE
  temp.dir = "output/DoAbsolute_mskcc_2020",
  results.dir = "output/DoAbsolute_mskcc_2020", nThread = 20,
  keepAllResult = TRUE, verbose = TRUE
)


# Take a check ------------------------------------------------------------


# Final output
df1 = fread("output/DoAbsolute_prad_su2c/DoAbsolute.wsx.ABSOLUTE.table.txt")
table(df1$`call status`)

df2 = fread("output/DoAbsolute_prad_su2c/summary/DoAbsolute.PP-calls_tab.txt")
table(df2$`call status`)

# Final output
df1 = fread("output/DoAbsolute_mskcc_2020/DoAbsolute.wsx.ABSOLUTE.table.txt")
table(df1$`call status`)

df2 = fread("output/DoAbsolute_mskcc_2020/summary/DoAbsolute.PP-calls_tab.txt")
table(df2$`call status`)



# Call copy number signatures ---------------------------------------------

library(sigminer)
library(tidyverse)
library(NMF)

# Set this per R session
options(sigminer.sex = "male", sigminer.copynumber.max = 20L)

cn_su2c = read_copynumber("output/DoAbsolute_prad_su2c/seg/")
cn_mskcc = read_copynumber("output/DoAbsolute_mskcc_2020/seg/")

save(cn_su2c, cn_mskcc, file = "output/CNV.validation.DoAbsolute.RData")

ncores <- 20

## Tally records
CNV.su2c.tally.W <- sig_tally(cn_su2c, method = "W", cores = ncores, feature_setting = CN.features)
CNV.mskcc.tally.W <- sig_tally(cn_mskcc, method = "W", cores = ncores, feature_setting = CN.features)


save(CNV.su2c.tally.W, file = "output/CNV.su2c.tally.W.RData")
save(CNV.mskcc.tally.W, file = "output/CNV.mskcc.tally.W.RData")


# ## Estimate signature number
# EST.su2c.W <- sig_estimate(CNV.su2c.tally.W$nmf_matrix,
#                                range = 2:8, nrun = 50, cores = ncores, use_random = TRUE,
#                                save_plots = FALSE, pConstant = 1e-9,
#                                verbose = TRUE
# )
# save(EST.su2c.W, file = "output/EST.su2c.W.RData")
#
# EST.mskcc.W <- sig_estimate(CNV.mskcc.tally.W$nmf_matrix,
#                            range = 2:8, nrun = 50, cores = ncores, use_random = TRUE,
#                            save_plots = FALSE, pConstant = 1e-9,
#                            verbose = TRUE
# )
# save(EST.mskcc.W, file = "output/EST.mskcc.W.RData")

## Extract signatures
## 5 signature, keep in line with previous analysis
load("output/CNV.su2c.tally.W.RData")
load("output/CNV.mskcc.tally.W.RData")

Sig.CNV.su2c = sig_extract(CNV.su2c.tally.W$nmf_matrix, n_sig = 5, nrun = 50,
                                                        cores = ncores, pConstant = 1e-9)
Sig.CNV.mskcc = sig_extract(CNV.mskcc.tally.W$nmf_matrix, n_sig = 5, nrun = 50,
                                                         cores = ncores, pConstant = 1e-9)
save(Sig.CNV.su2c, file = "output/Sig.CNV.su2c.RData")
save(Sig.CNV.mskcc, file = "output/Sig.CNV.mskcc.RData")


show_sig_profile(Sig.CNV.su2c, method = "W", normalize = "feature", style = "cosmic")
show_sig_profile(Sig.CNV.mskcc, method = "W", normalize = "feature", style = "cosmic")

show_sig_exposure(Sig.CNV.su2c, style = "cosmic", rm_space = TRUE)
show_sig_exposure(Sig.CNV.mskcc, style = "cosmic", rm_space = TRUE)

