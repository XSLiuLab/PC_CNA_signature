## Names used for FACETS
##
## The code is from Huimin and is not well described
## For reference only

## For sample name of dbGap studies
library("tidyverse")
sample <- mapping_df %>% select(gap_accession, subject_id, tumor_Run, normal_Run)
sample <- sample %>% separate(gap_accession, into = c("gap_accession", "accseion"), sep = 6)
sample <- sample %>% unite(subject_id, accseion, subject_id, sep = "-")
sample <- sample %>% select(subject_id, tumor_Run, normal_Run)


write.csv(sample, file = "C:/Users/lihm/Desktop/samplename", row.names = F, quote = F)



## For sample name of TCGA study
tcgapaired <- read.table("paired_sample.txt", sep = "\t", header = F, stringsAsFactors = F)
tcganame <- read.table("sampletcga", sep = ",", header = F, stringsAsFactors = F)


tumor_pool <- vector("character", length(12))
normal_pool <- vector("character", length(12))
j <- 1
k <- 1
for (i in tcganame$V1) {
  if (as.integer(substr(i, 19, 20)) < 10) {
    tumor_pool[j] <- c(i)
    j <- j + 1
  }
  else {
    normal_pool[k] <- c(i)
    k <- k + 1
  }
}


for (i in 1:length(tumor_pool)) {
  c <- unique(substr(tumor_pool, 6, 20))
  if (tcgapaired[i, 1] %in% c) {
    tcgapaired[i, 1] <- tumor_pool[which(c == tcgapaired[i, 1])]
  }
  i <- i + 1
}

for (i in 1:length(normal_pool)) {
  c <- unique(substr(normal_pool, 6, 20))
  if (tcgapaired[i, 2] %in% c) {
    tcgapaired[i, 2] <- normal_pool[which(c == tcgapaired[i, 2])]
  }
  i <- i + 1
}


write.csv(tcgapaired, file = "C:/Users/lihm/Desktop/tcgaid", row.names = F, quote = F)
