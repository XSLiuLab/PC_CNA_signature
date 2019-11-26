#PBS -N PBS_<sample>_seqz
#PBS -l nodes=1:ppn=1
#PBS -l walltime=70:00:00
#PBS -S /bin/bash
#PBS -j oe
#PBS -q normal_8

source activate R_36

out_dir=/public/home/wangshx/wangshx/PRAD_Sig
sseqz_dir=$out_dir"/small-seqz"
res_dir=$out_dir"/seqz_result"

Rscript $out_dir"/3-sequenza.R" <sample> $sseqz_dir/"<sample>.small.seqz.gz" $res_dir