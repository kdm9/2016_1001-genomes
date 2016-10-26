with open("metadata/srr_with_pheno_list.txt") as fh:
    RUNS = list(map(str.strip, fh.read().splitlines()))

onsuccess:
    shell("mail -s 'Workflow complete' kevin.murray@anu.edu.au < {log}")
onerror:
    shell("mail -s 'Workflow error' kevin.murray@anu.edu.au < {log}")


localrules: all, clean, sra

## BEGIN RULES
rule all:
    input:
        expand("data/reads/{run}.fastq.zst", run=RUNS),

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
        temp("data/tmp/reads/{run}.fastq.zst"),
    log:
        "data/log/dumpreads/{run}.log"
    shell:
        "(fastq-dump"
        "   --split-spot"
        "   --skip-technical"
        "   --stdout"
        "   --readids"
        "   --defline-seq '@$sn/$ri'"
        "   --defline-qual '+'"
        "   {input}"
        "| zstd -1 -o {output})"
        ">{log} 2>&1"

rule qcreads:
    input:
        "data/tmp/reads/{run}.fastq.zst",
    output:
        "data/reads/{run}.fastq.zst",
    log:
        "data/log/dumpreads/{run}.log"
    shell:
        "(( AdapterRemoval "
        "   --file1 <(zstdcat {input})"
        "   --output1 /dev/stdout"
        "   --interleaved"
        "   --combined-output"
        "   --settings /dev/null"
        "| sickle se"
        "   -t sanger"
        "   -q 28"
        "   -l 32"
        "   -g"
        "   -f /dev/stdin"
        "   -o >(zstd -19 -o {output}) )"
        "|| cp {input} {output} )"
        ">{log} 2>&1"
