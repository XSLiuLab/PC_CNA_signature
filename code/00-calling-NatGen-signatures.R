library(sigminer)
library(NMF)

data_maf = read_maf('/public/data/cbioPortal/signature/data_mutations_extended.txt')
maf.pre = sig_derive(data_maf, prefix = "chr", add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
maf.sig = sig_extract(maf.pre$nmf_matrix, n_sig = 5, nrun = 100, cores = 16)

get_sig_similarity(maf.sig)
sig.group = get_groups(Signature = maf.sig, method = "consensus")
table(sig.group$enrich_sig)
