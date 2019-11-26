source activate python3

ref_file=/public/home/liuxs/biodata/reference/genome/gdc_hg38/GRCh38.d1.vd1.fa
gc_file=/public/home/wangshx/wangshx/PRAD_Sig/hg38.gc50Base.wig.gz

mkdir -p $(dirname $gc_file)
# sequenza-utils version: 3.0.0
sequenza-utils gc_wiggle -w 50 --fasta $ref_file -o $gc_file
