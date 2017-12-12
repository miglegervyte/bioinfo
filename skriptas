#! /bin/sh 
# Skriptas sumappina visus meginius, esancius direktorijoj, kurie atitinka vardo sablona. 

bwa index chr2.fa


for NR in 1; ## tiek skaičių, kiek read'ų
do

readas=read${NR}.fastq

bwa aln chr2.fa $readas | bwa mem chr2.fa $readas | awk -F'\t' '$5 != "0"' - | samtools view -b -S - | bedtools bamtobed -i - > mapped-mapq.bed

done
