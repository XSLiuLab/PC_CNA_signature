#! /usr/bin/env Rscript

library("pctGCdata")
library("facets")
set.seed(1234)
rcmat = readSnpMatrix("/public/home/liuxs/ncbi/dbGaP-16533/copy/out/<sample>.out.gz")
xx = preProcSample(rcmat,gbuild = "hg38")
oo=procSample(xx,cval=150)
fit=emcncf(oo)

#plot
pdf("/public/home/liuxs/ncbi/dbGaP-16533/copy/facetdata_150/<sample>.pdf")
plotSample(x=oo,emfit=fit)
logRlogORspider(oo$out, oo$dipLogR)
while (!is.null(dev.list()))  dev.off()
#
save(fit,file = "/public/home/liuxs/ncbi/dbGaP-16533/copy/facetdata_150/<sample>.Rdata")
# output purity and ploidy -----
 purity=fit$purity
 purity=round(purity,2)
 ploidy=fit$ploidy
 ploidy=round(ploidy,1)
 output <- paste("<sample>", purity, ploidy, sep = "\t")
 write(output, "/public/home/liuxs/ncbi/dbGaP-16533/copy/facetdata_150/<sample>.txt", append = TRUE)
