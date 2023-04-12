#!/bin/bash

ls /mnt/project/Analysis_1/genotype/ > cram_run.txt
echo "IID   MeanCoverage" > tryptase_test.sample.tsv

for $cram in $(cat cram_run.txt)
do
    dx download Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/*/"${cram}".cram -o cram/
    dx download Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/*/"${cram}".cram.crai -o cram/

    samtools depth cram/${cram} | awk '{sum+=$3} END { print "${cram}   ",sum/NR}' >> tryptase_test.sample.tsv

    awk '{print "${cram}   ",$0}' >> tryptase_test.genotype.tsv
done