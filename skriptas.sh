#! /bin/sh 
# Tarkim, kad visi read'ai pavadinti 'read1.fastq', 'read2.fastq' ir pan. Referentinė seka - 'chr.fa'

bwa index chr.fa


for NR in 1; ## tiek skaičių, kiek read'ų
do

readas=read${NR}.fastq

bwa aln chr.fa $readas | bwa samse chr.fa - $readas | awk -F'\t' '$5 != "0"' - | samtools view -b -S - | samtools sort | bedtools bamtobed -i - > mapped-mapq.bed

done
