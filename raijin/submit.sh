logdir=raijin/log
mkdir -p $logdir

if [ -d ".kdmwrappers/.git" ]
then
    cd .kdmwrappers && git pull && cd ..
else
    git clone 'https://github.com/kdmurray91/snakemake-wrappers.git' .kdmwrappers
fi

QSUB="qsub -q {cluster.q} -l ncpus={threads} -l jobfs={cluster.jobfs}"
QSUB="$QSUB -l walltime={cluster.time} -l mem={cluster.mem}"
QSUB="$QSUB -l wd -o $logdir -e $logdir -P xe2"

snakemake                                \
    -j 300                               \
    --cluster-config raijin/cluster.yaml \
    --js raijin/jobscript.sh             \
    --rerun-incomplete                   \
    --keep-going                         \
    --cluster "$QSUB" $@
