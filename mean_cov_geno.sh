#!/bin/bash

ls /mnt/project/Analysis_1/genotype/ > cram_run.txt
echo "IID\tMeanCoverage" > tryptase_test.sample.tsv
samtools="/usr/local/bin/samtools"

for cram in $(cat cram_run.txt)
do
    dx download Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/10/"${cram}".cram -o cram/
    dx download Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/10/"${cram}".cram.crai -o cram/

    $samtools depth cram/${cram}.cram | awk -v cram_var="$cram" '{sum+=$3} END { print cram_var,"\t",sum/NR}' >> tryptase_test.sample.tsv

    #awk cram_var="$cram" '{print cram_var,"\t",$0}' /mnt/project/Analysis_1/genotype/${cram}/${cram}.output.txt >> tryptase_test.genotype.tsv
done