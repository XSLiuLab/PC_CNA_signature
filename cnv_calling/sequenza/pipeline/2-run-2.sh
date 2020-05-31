#!/bin/bash

pbs_dir='/public/home/wangshx/wangshx/PRAD_Sig/WES_calling_step1'

# 1. upload
loon upload WES_calling_step1/ $pbs_dir
# 2. submit
loon pbssub --remote --workdir $pbs_dir $pbs_dir 

# Deploy task directory to remote HPC
# The first argument must end with '/'
# loon pbsdeploy WES_calling_step1/ $pbs_dir -v --rsync