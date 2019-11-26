#!/bin/bash

pbs_dir='/public/home/wangshx/wangshx/PRAD_Sig/pbs'

# Test
loon upload deep_del_calling_pbs/*.pbs $pbs_dir -v
loon pbssub $pbs_dir \
    --remote -v \
    --workdir $pbs_dir
