export TERM=xterm
source /broad/software/scripts/useuse

use Java-1.7
use UGER
use Python-2.7
use Samtools

pip install -r requirements.txt --user --upgrade


cd /home/unix/weisburd/code/exac_readviz_scripts;

qsub -t 1-1000 \
-o /broad/hptmp/exac_readviz_backend/logs -e /broad/hptmp/exac_readviz_backend/logs \
-cwd -j y -V  -q short \
./run_python.sh python generate_HC_bams.py \
--process-variant-table  ` if [ "$#" -ne 0 ]; then echo --chrom $1; fi `
