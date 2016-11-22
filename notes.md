Data
====

- The 1135 Ath lines from SRA. https://www.ncbi.nlm.nih.gov/bioproject/273563
- Phenotype data from https://arapheno.1001genomes.org/


# 2016-11-18

Mitobim experiments

Using mitobim to try assembling the chloroplast of Col-0 (SRR1946065) and Cvi
(SRR1946067).

* Got chloroplast with `samtools faidx refs/tair10/TAIR10.fa ChrC >chloro.fa`

### MIRA mapping asm:

* Mira manifest.conf from mitobim tutorial in README on github

### Mitobim "quick":

* `mitobim -start 1 -end 3 -sample SRR1946065 -ref tair10 -readpool
  SRR1946065.fastq.gz --quick chloro.fa --verbose --pair
  --redirect_tmp $PBS_JOBFS`
* Assembled a full contig of the chloroplast

# 2016-11-22

Things seem to be crashing due to both walltime and RAM limits. Reducing
maximum number of iterations to 10 and increasing ram to 63G. Also setting the
priotity of Mitobim runs to -1, so other parts are run first. This should mean
we can analyse the chloroplast re-mapping pipeline for those who have
chloroplasts sooner.
