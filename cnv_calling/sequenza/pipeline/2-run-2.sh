#!/bin/bash

pbs_dir='/public/home/wangshx/wangshx/PRAD_Sig/pbs'

# Deploy task directory to remote HPC
# The first argument must end with '/'
loon pbsdeploy deep_del_calling_pbs/ $pbs_dir -v