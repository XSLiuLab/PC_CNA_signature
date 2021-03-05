# Single sample fitting testing for Macintyre method

## Preparation

1. Clone <https://bitbucket.org/britroc/cnsignatures/>
2. Set the clone repository as R project working directory
3. Check the single sample fitting with the repository's data and follow the clone repository README instruction

## Test

Obtain copy number data.

```R
library(tidyverse)
source("main_functions.R")
cn = readRDS("manuscript_Rmarkdown/data/britroc_absolute_copynumber.rds")
all_samples = colnames(cn)
```

Check with one sample.

```R
> cn2 = cn[, all_samples[1]]
> fts = extractCopynumberFeatures(cn2)
> mat = generateSampleByComponentMatrix(fts)
 Error in 1:ncol(posterior_sum) : argument of length 0 
8.
paste0(name, 1:ncol(posterior_sum)) at helper_functions.R#69
7.
calculateSumOfPosteriors(CN_features[["segsize"]], all_components[["segsize"]], 
    "segsize") at main_functions.R#210
6.
eval(quote(list(...)), env) 
5.
eval(quote(list(...)), env) 
4.
eval(quote(list(...)), env) 
3.
standardGeneric("cbind") 
2.
cbind(calculateSumOfPosteriors(CN_features[["segsize"]], all_components[["segsize"]], 
    "segsize"), calculateSumOfPosteriors(CN_features[["bp10MB"]], 
    all_components[["bp10MB"]], "bp10MB"), calculateSumOfPosteriors(CN_features[["osCN"]], 
    all_components[["osCN"]], "osCN"), calculateSumOfPosteriors(CN_features[["changepoint"]],  ... at main_functions.R#210
1.
generateSampleByComponentMatrix(fts) 
```

## Benchmark

To better understand how the code works, I write a function to return `1` if the analysis process is sucessess and `0` for failure.

```R
run_one = function(samples) {
  message("Samples: ", paste(samples, collapse = ", "))
  out = tryCatch({
    cn2 = cn[, samples]
    fts = extractCopynumberFeatures(cn2)
    mat = generateSampleByComponentMatrix(fts)
    expo = quantifySignatures(mat)
    1L
  }, error = function(e) {
    0L
  })
  return(out)
}
```

Then randomly select 1 to 10 samples from the data and run the analysis, repeat each analysis 100 times.

```R
library(furrr)
plan("multisession", workers = 10)
result = list()
for (i in 1:10) {
  message("Test sample number: ", i)
  result[[i]] = future_map_int(rerun(100, sample(all_samples, i)), run_one)
}
```

Check the result:

```R
> result
[[1]]
  [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
 [43] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
 [85] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

[[2]]
  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
 [43] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
 [85] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

[[3]]
  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
 [43] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
 [85] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

[[4]]
  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
 [43] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
 [85] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

[[5]]
  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
 [43] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
 [85] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

[[6]]
  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
 [43] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
 [85] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

[[7]]
  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
 [43] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
 [85] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

[[8]]
  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
 [43] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
 [85] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

[[9]]
  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
 [43] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
 [85] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

[[10]]
  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
 [43] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
 [85] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
```

## Conclusion

Based on data from original reference and the author's instruction, the Macintyre signature fitting analysis code works for >= 2 samples but not for single sample.
