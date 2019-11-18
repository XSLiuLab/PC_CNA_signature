# library(biomaRt)
#
# listEnsembl()
# ensembl <- useMart("ensembl")
# ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
#
# attributes = listAttributes(ensembl)
# attributes
#
# data = getBM(c("chromosome_name","start_position", "end_position", "ensembl_exon_id"), mart = ensembl)

# curl  -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/knownGene.txt.gz" | gunzip -c | \
# awk '{n=int($8); split($9,S,/,/);split($10,E,/,/); for(i=1;i<=n;++i) {printf("%s,%s,%s,%s,%s\n",$1,$2,$3,S[i],E[i]);} }'

system("curl ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz > data/hg38.gtf.gz")
gtf = data.table::fread("data/hg38.gtf.gz")
