#!/usr/bin/env bash

cd /home/unix/weisburd/code/exac_readviz_scripts/scripts/

d='/broad/hptmp/exac_readviz_backend2'


for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y;
do
	qsub -q long -o ${d}/logs/combine_chr${chrom}.log -cwd -j y -V ./run_combine_chrom.sh $chrom
done;
