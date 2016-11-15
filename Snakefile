REFERENCE='refs/tair10/TAIR10'

with open("metadata/srr_with_pheno_list.txt") as fh:
    RUNS = list(map(str.strip, fh.read().splitlines()))

# RUNS=['ERR1720545'] # testing tiny SRA

def kdmwrap(wrapper):
    local = ".kdmwrappers/" + wrapper
    if os.path.exists(local + "/wrapper.py"):
        return "file://" + local
    return "https://github.com/kdmurray91/snakemake-wrappers/raw/master/" + wrapper

localrules: all, clean, sra

## BEGIN RULES
rule all:
    input:
        expand("data/sketches/{run}.ct.gz", run=RUNS),
        expand("data/alignments/{run}.bam", run=RUNS),

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
        temp("data/tmp/reads/{run}.fastq"),
    log:
        "data/log/dumpreads/{run}.log"
    params:
        compressor="cat" # no compression
    wrapper:
        kdmwrap("sra/fastq-dump")

rule qcreads:
    input:
        "data/tmp/reads/{run}.fastq",
    output:
        "data/reads/{run}.fastq.gz",
    log:
        "data/log/qcreads/{run}.log"
    params:
        extra='-q 28 -l 32 -z',
        qual_type="sanger",
    wrapper:
        kdmwrap("sickle/se")

rule bwamem:
    input:
        ref=REFERENCE,
        sample=["data/reads/{run}.fastq.gz",],
    output:
        temp("data/tmp/aln_unsort/{run}.bam"),
    log:
        "data/logs/align/{run}.log"
    params:
        rgid=lambda wc: wc.run,
        extra="-p",
    threads: 16
    wrapper:
        kdmwrap("bwa/mem")


rule bamsort:
    input:
        "data/tmp/aln_unsort/{run}.bam"
    output:
        "data/alignments/{run}.bam",
        "data/alignments/{run}.bam.bai",
    log:
        "data/logs/bamsort/{run}.log"
    params:
        mem='3G',
        tmpdir='${TMPDIR:-/tmp}'
    threads: 8
    wrapper:
        kdmwrap("samtools/sortidx")


rule sketchreads:
    input:
        reads="data/reads/{run}.fastq.gz",
    output:
        ct="data/sketches/{run}.ct.gz",
        tsv="data/sketches/{run}.ct.gz.info.tsv",
        inf="data/sketches/{run}.ct.gz.info",
    log:
        "data/log/sketchreads/{run}.log"
    threads:
        8
    params:
        mem=12e9,
        k=31,
    shell:
        'load-into-counting.py'
        '    -M {params.mem}'
        '    -k {params.k}'
        '    -s tsv'
        '    -b'
        '    -f'
        '    -T {threads}'
        '    {output.ct}'
        '    {input}'
        '>{log} 2>&1'
