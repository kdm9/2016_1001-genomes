logdir=raijin/log
mkdir -p $logdir
QSUB="qsub -q {cluster.q} -l ncpus={threads} -l jobfs={cluster.jobfs}"
QSUB="$QSUB -l walltime={cluster.time} -l mem={cluster.mem}"
QSUB="$QSUB -l wd -o $logdir -e $logdir -P xe2"

snakemake                                \
    --keep-going                         \
    -j 250                               \
    --cluster-config raijin/cluster.yaml \
    --js raijin/jobscript.sh             \
    --cluster "$QSUB"
