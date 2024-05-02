#!/bin/bash

#echo -e "IID\tMeanCoverage" > tryptase_${1}.sample.tsv
input="bulk"${1}"cram_new.txt"

while read path
do
        cram_name=$(basename "$path" .cram)
        
        cat tsv/${cram_name}.sample.tsv >> tryptase_${1}.sample.tsv
        cat tsv/${cram_name}.genotype.tsv >> tryptase_${1}.genotype.tsv

done < $input

dx upload tryptase_${1}.sample.tsv --path /Analysis_1/tsv/tryptase_${1}.sample.tsv
dx upload tryptase_${1}.genotype.tsv --path /Analysis_1/tsv/tryptase_${1}.sample.tsv
