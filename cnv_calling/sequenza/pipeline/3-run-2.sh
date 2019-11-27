#!/bin/bash

pbs_dir='/public/home/wangshx/wangshx/PRAD_Sig/pbs_step2'

# upload R script
loon upload 3-sequenza.R '/public/home/wangshx/wangshx/PRAD_Sig'
# Deploy task directory to remote HPC
# The first argument must end with '/'
loon pbsdeploy deep_del_calling_pbs_step2/ $pbs_dir -v