#!/bin/bash

pbs_dir='/public/home/wangshx/wangshx/PRAD_Sig/WES_calling_step2'

# upload R script
loon upload 3-sequenza.R '/public/home/wangshx/wangshx/PRAD_Sig'

# 1. upload
loon upload WES_calling_step2/ $pbs_dir
# 2. submit
loon pbssub --remote --workdir $pbs_dir $pbs_dir 