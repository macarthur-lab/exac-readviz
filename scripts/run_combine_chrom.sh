#!/usr/bin/env bash

# runs combine on a single chrom

source /broad/software/scripts/useuse
use Python-2.7
use Java-1.7

i=$1

d='/broad/hptmp/exac_readviz_backend2'

cd /home/unix/weisburd/code/exac_readviz_scripts

echo Processing chrom $i
python2.7 -u -m scripts.combine_bams --chrom $i -k1   0 -k2 249 2>&1 | grep -v Warning &    # 2>&1 > ${d}/logs/chr${i}_000.log &
python2.7 -u -m scripts.combine_bams --chrom $i -k1 250 -k2 499 2>&1 | grep -v Warning &   # 2>&1 > ${d}/logs/chr${i}_250.log &
python2.7 -u -m scripts.combine_bams --chrom $i -k1 500 -k2 749 2>&1 | grep -v Warning &   # 2>&1 > ${d}/logs/chr${i}_500.log &
python2.7 -u -m scripts.combine_bams --chrom $i -k1 750 -k2 999 2>&1 | grep -v Warning &   # 2>&1 > ${d}/logs/chr${i}_750.log &
echo Launched 4 processes

for job in `jobs -p`
do
    echo Waiting for $job
    wait $job || let "FAIL+=1"
done

