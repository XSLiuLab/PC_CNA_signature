load(file = "output/Sig.PRAD_TCGA_plus_dbGap_rm_hyper.RData")
# Method W
load(file = "output/Sig.CNV.seqz.W.5.RData")

load(file = "output/Sig.CNV.seqz.W.RData")
load(file = "output/Sig.CNV.facets.W.RData")
# Method M
load(file = "output/Sig.CNV.facets.M.RData")
load(file = "output/Sig.CNV.seqz.M.RData")
load(file = "output/Sig.CNV.facets.M.ref.seqz.RData")

# Derive data
load(file = "output/CNV.seqz.derive.W.RData")
load(file = "output/CNV.facets.derive.W.RData")
load(file = "output/CNV.seqz.derive.M.RData")
load(file = "output/CNV.facets.derive.M.RData")
load(file = "output/CNV.facets.derive.M.ref.seqz.RData")


# Signature profile -------------------------------------------------------

show_sig_profile(Sig.CNV.seqz.W.5, method = "W", normalize = "feature", x_label_angle = 90)
show_sig_profile(Sig.CNV.seqz.W, method = "W", normalize = "feature", x_label_angle = 90)
show_sig_profile(Sig.CNV.facets.W, method = "W", normalize = "feature", x_label_angle = 90)


show_sig_profile(Sig.CNV.seqz.M,
                 method = "M", normalize = "feature",
                 params = CNV.seqz.derive.M$parameters, y_expand = 1.5,
                 set_gradient_color = FALSE,
                 x_label_angle = 90
)

show_sig_profile(Sig.CNV.facets.M,
                 method = "M", normalize = "feature",
                 params = CNV.facets.derive.M$parameters, y_expand = 1.5,
                 set_gradient_color = FALSE,
                 x_label_angle = 90
)


show_sig_profile(Sig.CNV.facets.M.ref.seqz,
                 method = "M", normalize = "feature",
                 params = CNV.facets.derive.M.ref.seqz$parameters, y_expand = 1.5,
                 set_gradient_color = FALSE,
                 x_label_angle = 90
)


show_sig_profile(Sig.SNV, mode = "mutation")

# Signature exposure ------------------------------------------------------
show_sig_exposure(Sig.SNV, rm_space = T, cutoff = 2000)

show_sig_exposure(Sig.CNV.seqz.W, rm_space = T)
show_sig_exposure(Sig.CNV.facets.W, rm_space = T)


dd <- get_groups(Sig.CNV.seqz.W, method = "consensus")
table(dd$enrich_sig)

dd <- get_groups(Sig.CNV.seqz.W.5, method = "consensus")
table(dd$enrich_sig)

# Signature similarity ----------------------------------------------------

get_sig_similarity(Sig.CNV.seqz.W, Sig.CNV.facets.W, normalize = "feature")

get_sig_similarity(Sig.CNV.seqz.W.5, Sig.CNV.facets.W, normalize = "feature")
get_sig_similarity(Sig.CNV.seqz.W.5, Sig.CNV.seqz.W, normalize = "feature")

get_sig_similarity(Sig.CNV.seqz.M, Sig.CNV.facets.M.ref.seqz, normalize = "feature")

get_sig_similarity(Sig.SNV)
get_sig_similarity(Sig.SNV, sig_db = "SBS")



# Check signature number survey -------------------------------------------

load("output/EST.seqz.W.all.RData")
load("output/EST.facets.W.all.RData")
load("output/EST.PRAD_TCGA_plus_dbGap_Maf_rm_hyper.RData")

show_sig_number_survey(EST.seqz.W.all)
show_sig_number_survey2(EST.seqz.W.all$survey)

show_sig_number_survey(EST.facets.W.all)
show_sig_number_survey2(EST.facets.W.all$survey)

show_sig_number_survey(EST.Maf.rm_hyper)
show_sig_number_survey2(EST.Maf.rm_hyper$survey, EST.Maf.rm_hyper$survey.random)
