#!/bin/bash

pbs_dir='/public/home/wangshx/wangshx/PRAD_Sig/pbs_step2'

# Test
loon upload 3-sequenza.R '/public/home/wangshx/wangshx/PRAD_Sig/'
loon upload deep_del_calling_pbs_step2/*.pbs $pbs_dir -v
loon pbssub $pbs_dir \
    --remote -v \
    --workdir $pbs_dir