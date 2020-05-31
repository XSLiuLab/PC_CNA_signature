#! /bin/bash
# Generate batch pbs

loon pbsgen -t 3-sequenza-1.sh -s WES_samples.csv \
    -m 3-seqz-wes.mapfile -o WES_calling_step2