if [ ! -d /humgen/atgu1/fs03 ]; then
    echo "/humgen/atgu1/fs03 not found. Exiting.."
    exit -1
fi

export TERM=xterm
source /broad/software/scripts/useuse
reuse UGER
reuse Python-3.4
reuse Python-2.7

echo "run_python.sh: user=$USER"
echo "run_python.sh: host=$HOST"
echo "run_python.sh: load="`uptime`

$@

