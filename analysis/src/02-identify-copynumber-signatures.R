
# Load packages -----------------------------------------------------------

library(sigminer)
library(tidyverse)
library(NMF)

# Set this per R session
options(sigminer.sex = "male", sigminer.copynumber.max = 20L)

# Generate CopyNumber object ----------------------------------------------
# Note: we used two name system for FACETS and Sequenza data

CNV.facets <- read_copynumber("data/CNV_from_dbGAP_PLUS_TCGA_WES_CVAL150.tsv",
  genome_build = "hg38",
  complement = FALSE, verbose = TRUE
)
# find 141-10 only have segments in chr1
CNV.facets <- subset(CNV.facets, subset = !sample %in% "141-10")

save(CNV.facets, file = "output/CNV.facets.RData")

CNV.seqz <- read_copynumber("data/CNV_from_sequenza.tsv",
  genome_build = "hg38",
  complement = FALSE, verbose = TRUE
)

# remove WCMC160-SRR3146971 with only one CNV
CNV.seqz <- subset(CNV.seqz, subset = !sample %in% "WCMC160-SRR3146971")
save(CNV.seqz, file = "output/CNV.seqz.RData")

load(file = "output/CNV.seqz.RData")
load(file = "output/CNV.facets.RData")
# tally copy number features ---------------------------------------------

ncores <- 20

##
## W method
##

# Use classfication method devised by me ("W")
#  85.984 s
CNV.seqz.tally.W <- sig_tally(CNV.seqz, method = "W", cores = ncores, feature_setting = CN.features)
save(CNV.seqz.tally.W, file = "output/CNV.seqz.tally.W.RData")

CNV.facets.tally.W <- sig_tally(CNV.facets, method = "W", cores = ncores, feature_setting = CN.features)
save(CNV.facets.tally.W, file = "output/CNV.facets.tally.W.RData")

##
## M method
##

# Use classfication method from Macintyre et al ("M")
system.time(
  CNV.seqz.tally.M <- sig_tally(CNV.seqz, method = "M", cores = ncores, nrep = 3)
)
# 5126.994s
save(CNV.seqz.tally.M, file = "output/CNV.seqz.tally.M.RData")

system.time(
  CNV.facets.tally.M <- sig_tally(CNV.facets, method = "M", cores = ncores, nrep = 3)
)
save(CNV.facets.tally.M, file = "output/CNV.facets.tally.M.RData")

# Use components from sequenza as reference
CNV.facets.tally.M.ref.seqz <- sig_tally(CNV.facets,
  method = "M",
  reference_components = CNV.seqz.tally.M$components,
  cores = ncores
)

save(CNV.facets.tally.M.ref.seqz, file = "output/CNV.facets.tally.M.ref.seqz.RData")



# Estimate number of copy number signatures -------------------------------
ncores <- 20

#
EST.seqz.W.all <- sig_estimate(CNV.seqz.tally.W$nmf_matrix,
  range = 2:12, nrun = 50, cores = ncores, use_random = TRUE,
  save_plots = FALSE,
  verbose = TRUE
)
save(EST.seqz.W.all, file = "output/EST.seqz.W.all.RData")

#
EST.facets.W.all <- sig_estimate(CNV.facets.tally.W$nmf_matrix,
  range = 2:12, nrun = 50, cores = ncores, use_random = TRUE,
  save_plots = FALSE, pConstant = 1e-9,
  verbose = TRUE
)
save(EST.facets.W.all, file = "output/EST.facets.W.all.RData")

#
EST.seqz.M <- sig_estimate(CNV.seqz.tally.M$nmf_matrix,
  range = 2:12, nrun = 50, cores = ncores, use_random = TRUE,
  save_plots = FALSE,
  verbose = TRUE
)
save(EST.seqz.M, file = "output/EST.seqz.M.RData")

#
EST.facets.M <- sig_estimate(CNV.facets.tally.M$nmf_matrix,
  range = 2:12, nrun = 50, cores = ncores, use_random = TRUE,
  save_plots = FALSE,
  verbose = TRUE
)
save(EST.facets.M, file = "output/EST.facets.M.RData")

#
EST.facets.M.ref.seqz <- sig_estimate(CNV.facets.tally.M.ref.seqz$nmf_matrix,
  range = 2:12, nrun = 50, cores = ncores, use_random = TRUE,
  save_plots = FALSE,
  verbose = TRUE
)
save(EST.facets.M.ref.seqz, file = "output/EST.facets.M.ref.seqz.RData")


# Check the proper signature number -----------------------------------------

load("output/EST.seqz.W.all.RData")
load("output/EST.facets.W.all.RData")
load("output/EST.seqz.M.RData")
load("output/EST.facets.M.RData")
load("output/EST.facets.M.ref.seqz.RData")

# show_sig_number_survey(EST.seqz.W)
# show_sig_number_survey(EST.facets.W)

show_sig_number_survey(EST.seqz.W.all)
EST.seqz.W.all$survey_plot
show_sig_number_survey(EST.facets.W.all)
EST.facets.W.all$survey_plot

show_sig_number_survey(EST.seqz.M)
show_sig_number_survey(EST.facets.M)

show_sig_number_survey(EST.facets.M.ref.seqz)

# Extract copy number signatures ------------------------------------------

load(file = "output/CNV.seqz.tally.W.RData")
load(file = "output/CNV.facets.tally.W.RData")
load(file = "output/CNV.seqz.tally.M.RData")
load(file = "output/CNV.facets.tally.M.RData")
load(file = "output/CNV.facets.tally.M.ref.seqz.RData")


Sig.CNV.seqz.W <- sig_extract(CNV.seqz.tally.W$nmf_matrix, n_sig = 5, nrun = 50, cores = ncores)
## Keep in line with 5 signatures
Sig.CNV.facets.W <- sig_extract(CNV.facets.tally.W$nmf_matrix, n_sig = 5, nrun = 50, cores = ncores, pConstant = 1e-9)
Sig.CNV.seqz.M <- sig_extract(CNV.seqz.tally.M$nmf_matrix, n_sig = 5, nrun = 50, cores = ncores)
Sig.CNV.facets.M <- sig_extract(CNV.facets.tally.M$nmf_matrix, n_sig = 5, nrun = 50, cores = ncores)
Sig.CNV.facets.M.ref.seqz <- sig_extract(CNV.facets.tally.M.ref.seqz$nmf_matrix, n_sig = 5, nrun = 50, cores = ncores)

save(Sig.CNV.seqz.W, file = "output/Sig.CNV.seqz.W.RData")
save(Sig.CNV.facets.W, file = "output/Sig.CNV.facets.W.RData")
save(Sig.CNV.seqz.M, file = "output/Sig.CNV.seqz.M.RData")
save(Sig.CNV.facets.M, file = "output/Sig.CNV.facets.M.RData")
save(Sig.CNV.facets.M.ref.seqz, file = "output/Sig.CNV.facets.M.ref.seqz.RData")


# Some checks and analysis ------------------------------------------------

get_sig_similarity(Sig.CNV.seqz.W, Sig.CNV.facets.W, normalize = "feature")

show_sig_profile(Sig.CNV.seqz.W, method = "W", normalize = "feature", style = "cosmic")
show_sig_profile(Sig.CNV.seqz.W, method = "W", normalize = "feature", style = "cosmic",
                 filters = "BoChr", palette = use_color_style("cosmic")[8])
show_sig_profile(Sig.CNV.facets.W, method = "W", normalize = "feature", style = "cosmic")

show_sig_exposure(Sig.CNV.seqz.W, rm_space = T)
show_sig_exposure(Sig.CNV.facets.W, rm_space = T)

# show_sig_profile(Sig.CNV.seqz.W, method = "W",
#                  normalize = "feature", style = "cosmic",
#                  sig_orders = paste0("Sig",
#                                      sigminer:::helper_sort_signature(Sig.CNV.seqz.W$Signature.norm)))


## Test NC50

cn_list = sigminer:::get_cnlist(CNV.seqz)
nc50_dt = sigminer:::getNC50(cn_list)
hist(nc50_dt$value)
table(nc50_dt$value)
