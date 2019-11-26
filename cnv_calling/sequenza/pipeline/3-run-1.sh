#! /bin/bash
# Generate batch pbs

loon pbsgen -t 3-sequenza-1.sh -s deep_del_samples.csv \
    -m 3-seqz-wes.mapfile -o deep_del_calling_pbs_step2