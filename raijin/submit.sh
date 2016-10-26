logdir=data/log/cluster
mkdir -p $logdir
qsubargs="-l wd -o $logdir -e $logdir -P xe2"

snakemake                                \
    --keep-going                         \
    -j 390                               \
    --cluster-config raijin/cluster.json \
    --js raijin/jobscript.sh             \
    --cluster "qsub -q {cluster.q} -l ncpus={threads} -l walltime={cluster.time} -l mem={cluster.mem} $qsubargs"

