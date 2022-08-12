#!/bin/bash
set -x
echo project: $1
echo wd: $2
echo uuid_list: $3


qsub -q ngs4G \
     -P MST108173 \
	 -W group_list=MST108173 \
	 -N combine_count \
	 -l select=1:ncpus=1 \
	 -l place=pack \
	 -o $1 \
	 -e $1 \
	 -- /project/GP1/chialang1220/kpc/combine_counts.sh "$1" "$2" "$3"
