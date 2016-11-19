REFERENCES={
    "tair10": lambda wc: 'refs/tair10/TAIR10.fa',
    "chloro": lambda wc: "data/chloroplasts/{}.fasta".format(wc.run)
}

with open("metadata/srr_with_pheno_list.txt") as fh:
    RUNS = list(map(str.strip, fh.read().splitlines()))

# RUNS=['ERR1720545'] # testing tiny SRA
# RUNS=['SRR1946065', 'SRR1946067'] # Col and CVI resp.

def kdmwrap(wrapper):
    local = ".kdmwrappers/" + wrapper
    if os.path.exists(local + "/wrapper.py"):
        return "file://" + local
    return "https://github.com/kdmurray91/snakemake-wrappers/raw/master/" + wrapper

shell.executable("/bin/bash")
# bash "safe mode"
shell.prefix("set -euo pipefail; ")

localrules: all, clean, sra

## BEGIN RULES
rule all:
    input:
        expand("data/sketches/{run}.ct.gz", run=RUNS),
        expand("data/alignments/tair10/{run}.bam", run=RUNS),
        expand("data/alignments/chloro/{run}.bam", run=RUNS),
        expand("data/chloroplasts/{run}.fasta", run=RUNS),

rule clean:
    shell:
        "rm -rf data .snakemake"

rule sra:
    output:
        "data/sra/{run}.sra",
    log:
        "data/logs/getrun/{run}.log"
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
        "data/logs/dumpreads/{run}.log"
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
        "data/logs/qcreads/{run}.log"
    params:
        extra='-q 28 -l 32 -z',
        qual_type="sanger",
    wrapper:
        kdmwrap("sickle/se")


rule sketchreads:
    input:
        reads="data/reads/{run}.fastq.gz",
    output:
        ct="data/sketches/{run}.ct.gz",
        tsv="data/sketches/{run}.ct.gz.info.tsv",
        inf="data/sketches/{run}.ct.gz.info",
    log:
        "data/logs/sketchreads/{run}.log"
    threads:
        4
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

rule mitobim:
    input:
        reads="data/reads/{run}.fastq.gz",
        bait="refs/tair10/chloroplast.fasta",
    output:
        tarball="data/chloroplasts/{run}.tar.gz",
        fasta="data/chloroplasts/{run}.fasta"
    log:
        "data/logs/mitobim/{run}.log"
    threads:
        1
    shell:
        "(export WKDIR=$PWD;"
        " rm -rf $PBS_JOBFS/{wildcards.run}"
        " && mkdir $PBS_JOBFS/{wildcards.run}"
        " && cd $PBS_JOBFS/{wildcards.run}"
        " && mitobim"
        "   -start 1"
        "   -end 10"
        "   -sample {wildcards.run}"
        "   -ref bait"
        "   -readpool $WKDIR/{input.reads}"
        "   --quick $WKDIR/{input.bait}"
        "   --verbose"
        "   --pair"
        "   --symlink-final"
        " && cp final_iteration/*_assembly/*_d_results/*_out_AllStrains.unpadded.fasta"
        "      $WKDIR/{output.fasta}"
        " && (find -type d -and -name \*_chkpt -or -name \*_tmp | xargs rm -rf )"
        " && tar cvzf $WKDIR/{output.tarball} ."
        "; unset WKDIR) >{log} 2>&1"

rule bwaidx:
    input:
        "{fasta,.+(fasta|fa)}"
    output:
        "{fasta}.bwt"
    threads: 1
    shell:
        "bwa index {input} && "
        "samtools faidx {input}"


rule bwamem:
    input:
        ref=lambda w: REFERENCES[w.ref](w),
        sample=["data/reads/{run}.fastq.gz",],
        refidx=lambda w: REFERENCES[w.ref](w) + '.bwt',
    output:
        temp("data/tmp/aln_unsort/{ref}/{run}.bam"),
    log:
        "data/logs/align/{ref}/{run}.log"
    params:
        rgid=lambda wc: wc.run,
        extra="-p",
    threads: 16
    wrapper:
        kdmwrap("bwa/mem")


rule bamsort:
    input:
        "data/tmp/aln_unsort/{ref}/{run}.bam"
    output:
        "data/alignments/{ref}/{run}.bam",
        "data/alignments/{ref}/{run}.bam.bai",
    log:
        "data/logs/bamsort/{ref}/{run}.log"
    params:
        mem='20G',
        tmpdir='${TMPDIR:-/tmp}'
    threads: 1
    wrapper:
        kdmwrap("samtools/sortidx")
