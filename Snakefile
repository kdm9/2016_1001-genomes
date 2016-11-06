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
        expand("data/sketches/{run}.ct", run=RUNS),

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
        "| zstd --quiet -1 -o {output})"
        ">{log} 2>&1"

rule qcreads:
    input:
        "data/tmp/reads/{run}.fastq.zst",
    output:
        "data/reads/{run}.fastq.zst",
    log:
        "data/log/qcreads/{run}.log"
    threads:
        4
    shell:
        "(( AdapterRemoval "
        "   --file1 <(zstdcat --quiet {input})"
        "   --output1 /dev/stdout"
        "   --interleaved"
        "   --combined-output"
        "   --settings /dev/null"
        "| sickle se"
        "   -t sanger"
        "   -q 28"
        "   -l 32"
        "   -f /dev/stdin"
        "   -o >(pzstd -p {threads} -f -15 --quiet -o {output}) )"
        "|| (zstdcat --quiet {input} | pzstd -f -p {threads} -15 --quiet -o {output}))"
        ">{log} 2>&1"

rule sketchreads:
    input:
        "data/reads/{run}.fastq.zst",
    output:
        "data/sketches/{run}.ct",
    log:
        "data/log/sketchreads/{run}.log"
    threads:
        4
    params:
        mem=12e9,
        k=31,
    shell:
        '(zstdcat --quiet {input} '
        '| load-into-counting.py '
        '    -M {params.mem} '
        '    -k {params.k} '
        '    -s tsv '
        '    -b '
        '    -f '
        '    -T {threads} '
        '    {output} '
        '    - '
        ') >{log} 2>&1 '
