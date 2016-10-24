with open("metadata/srr_with_pheno_list.txt") as fh:
    RUNS = list(map(str.strip, fh.read().splitlines()))

onsuccess:
    shell("mail -s 'Workflow complete' kevin.murray@anu.edu.au < {log}")
onerror:
    shell("mail -s 'Workflow error' kevin.murray@anu.edu.au < {log}")


## BEGIN RULES
rule all:
    input:
        expand("data/reads/{run}_il.fastq.gz", run=RUNS),

rule clean:
    shell:
        "rm -rf data .snakemake"

rule sra:
    output:
        "data/sra/{run}.sra",
    log:
        "data/log/getrun/{run}.log"
    shell:
        "get-run.py"
        "   -d data/sra/"
        "   -s "
        "   -i {wildcards.run}"
        " >{log} 2>&1"

rule dumpreads:
    input:
        "data/sra/{run}.sra",
    output:
        "data/reads/{run}_il.fastq.gz",
    log:
        "data/log/dumpreads/{run}.log"
    shell:
        "( fastq-dump"
        "   --split-spot"
        "   --skip-technical"
        "   --stdout"
        "   --readids"
        "   --defline-seq '@$sn/$ri'"
        "   --defline-qual '+'"
        "   {input}"
        "| AdaptorRemoval "
        "   --file1 /dev/stdin"
        "   --output1 /dev/stdout"
        "   --interleaved"
        "   --combined-output"
        "| sickle se"
        "   -t sanger"
        "   -q 28"
        "   -l 32"
        "   -f /dev/stdin"
        "   -o /dev/stdout"
        "| gzip >{output})"
        ">{log} 2>&1"
