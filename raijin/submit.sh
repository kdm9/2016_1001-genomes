mkdir -p data/log/cluster

snakemake                                \
    -j 150                               \
    --cluster-config raijin/cluster.json \
    --js raijin/jobscript.sh             \
    --cluster 'qsub -q {cluster.queue} -l ncpus={threads} -l walltime={cluster.time} -l mem={cluster.mem} {cluster.miscargs}'

