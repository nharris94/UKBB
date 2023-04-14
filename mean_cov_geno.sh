#!/bin/bash

ls /mnt/project/Analysis_1/genotype/ > cram_run.txt
echo -e "IID\tMeanCoverage" > tryptase_test.sample.tsv
samtools="/usr/local/bin/samtools"
i=0

for cram in $(cat cram_run.txt)
do
    dx download Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/10/"${cram}".cram -o cram/
    dx download Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/10/"${cram}".cram.crai -o cram/

    $samtools depth cram/${cram}.cram | awk -v cram_var="$cram" 'OFS="\t" {sum+=$3} END { print cram_var,sum/NR}' >> tryptase_test.sample.tsv

    awk -v cram_var="$cram" 'OFS="\t" {print cram_var,$0}' /mnt/project/Analysis_1/genotype/${cram}/${cram}.output.txt >> tryptase_test.genotype.tsv
    
    rm cram/${cram}.*
    
    #after running certain number of jobs, terminate
    #((i=i+1))
    #if [[ $i == 5 ]]
    #then
    #    break
    #fi
done