#!/bin/bash
# Download small-seqz files with error in sequenza fit

work_dir=/public/home/wangshx/wangshx/PRAD_Sig
result_dir=$work_dir"/small-seqz"
store_dir=/home/wsx/projects/prad_signature/cnv_calling/sequenza/seqz_error_input

files=$(loon run cat $work_dir"/wes_error_input_files.txt")

for fl in $files
do
    echo $fl
    loon download $result_dir"/"$fl $store_dir #--dry
done