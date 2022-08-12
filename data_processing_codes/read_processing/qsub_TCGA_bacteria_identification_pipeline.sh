#!/bin/bash
set -x
#echo wd: $1
#echo uuid_list: $2

for uuid in $(xargs < $2)
do
	qsub -q ngs8G \
	     -P MST108173 \
		 -W group_list=MST108173 \
		 -N sRNAnalyzer_1_sample_test \
		 -l select=1:ncpus=2 \
		 -l place=pack \
		 -o $1 \
		 -e $1 \
		 -- /project/GP1/chialang1220/kpc/TCGA_bacteria_identification_pipeline.sh "$1" "${uuid}"
done
