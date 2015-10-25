export TERM=xterm
source /broad/software/scripts/useuse
reuse UGER
reuse Python-3.4

echo "run_python.sh: user=$USER"
echo "run_python.sh: host=$HOST"
echo "run_python.sh: load="`uptime`

$@

