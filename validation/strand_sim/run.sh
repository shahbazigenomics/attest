#!/usr/bin/env bash
# Real STAR and featureCounts on the simulated libraries.
set -euo pipefail
cd "$(dirname "$0")"
python3 simulate.py
mkdir -p idx star fc
STAR --runMode genomeGenerate --genomeDir idx --genomeFastaFiles ref/genome.fa \
     --sjdbGTFfile ref/genes.gtf --sjdbOverhang 74 --genomeSAindexNbases 9 \
     --runThreadN 4 > /dev/null
for fq in fastq/*.fq.gz; do
  n=$(basename "$fq" .fq.gz)
  STAR --genomeDir idx --readFilesIn "$fq" --readFilesCommand zcat \
       --outSAMtype BAM Unsorted --quantMode GeneCounts --runThreadN 4 \
       --outFileNamePrefix "star/$n." > /dev/null
done
# featureCounts on every library with each strand setting, as a user might
bams=$(ls star/*.Aligned.out.bam)
for s in 0 1 2; do
  featureCounts -s $s -a ref/genes.gtf -o fc/counts_s$s.txt -T 4 $bams > /dev/null 2>&1
done
ls star/*ReadsPerGene.out.tab fc/*.summary
# one featureCounts run per library type, as the package's example files
for s in 0 1 2; do featureCounts -s $s -a ref/genes.gtf -o fc/rev_s$s.txt -T 4 star/reverse_s*.Aligned.out.bam > /dev/null 2>&1; done
featureCounts -s 2 -a ref/genes.gtf -o fc/unstr_s2.txt  -T 4 star/unstranded_s*.Aligned.out.bam > /dev/null 2>&1
featureCounts -s 2 -a ref/genes.gtf -o fc/revbg_s2.txt  -T 4 star/reverse_bg_s*.Aligned.out.bam > /dev/null 2>&1
