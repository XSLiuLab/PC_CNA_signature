#! /bin/bash
# Run under /public/home/wangshx/wangshx/PRAD_Sig

pbs_dir=WES_calling_step2
res_dir=../seqz_wes_result

cd $pbs_dir
for f in $(ls *.pbs | xargs -n1)
do
	echo "=> Checking file $f..."
	name=$(echo $f | sed 's/\.pbs//')
	check=$(ls -l $res_dir | grep $name)
	if [ -z "$check" ]
	then 
		echo $name".small_filter.seqz.gz" >> ../wes_error_input_files.txt
	else
		echo "==> Result has been generated for $f"
		echo "==> Skipping..."
	fi
done