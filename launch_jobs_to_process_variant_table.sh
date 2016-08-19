export TERM=xterm
source /broad/software/scripts/useuse

use Java-1.7
use UGER
use Python-2.7
use Samtools

echo "Launching array with $1 tasks";

pip install --user --upgrade pip wheel

cd /home/unix/weisburd/code/exac_readviz_scripts && \
pip install -r requirements.txt --user --upgrade  && \
python2.7 -u parallelize.py -isize 1000000 -n $1 -L /home/unix/weisburd/code/exac_readviz_scripts/scripts/data/25_hc_runs.intervals python2.7 -u -m pipeline.compute_HC_bams_from_sample_table

# /seq/references/Homo_sapiens_assembly19/v1/variant_calling/exome_calling_regions.v1.interval_list
#qsub -t 1-1000 -o /broad/hptmp/exac_readviz_backend/logs -e /broad/hptmp/exac_readviz_backend/logs -cwd -j y -V  -q short \
#./run_python.sh python generate_HC_bams.py --process-variant-table  ` if [ "$#" -ne 0 ]; then echo --chrom $1; fi `

