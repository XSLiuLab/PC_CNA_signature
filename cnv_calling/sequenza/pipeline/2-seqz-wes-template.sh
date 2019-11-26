#PBS -N PBS_<sample>_seqz
#PBS -l nodes=1:ppn=1
#PBS -l walltime=70:00:00
#PBS -S /bin/bash
#PBS -j oe
#PBS -q normal_8

source activate python3

ref_file=/public/home/liuxs/biodata/reference/genome/gdc_hg38/GRCh38.d1.vd1.fa
gc_file=/public/home/wangshx/wangshx/PRAD_Sig/hg38.gc50Base.wig.gz

out_dir=/public/home/wangshx/wangshx/PRAD_Sig
seqz_dir=$out_dir"/seqz"
sseqz_dir=$out_dir"/small-seqz"

mkdir -p $seqz_dir
mkdir -p $sseqz_dir

#tumor=<tumor>
#normal=<normal>

sequenza-utils bam2seqz -n <normal> -t <tumor> \
	--fasta $ref_file -gc $gc_file \
	-o $seqz_dir/"<sample>.seqz.gz"

sequenza-utils seqz_binning --seqz $seqz_dir/"<sample>.seqz.gz" \
	-w 50 -o $sseqz_dir/"<sample>.small.seqz.gz"


