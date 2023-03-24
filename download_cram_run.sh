#!/bin/bash

i=0
input="bulk"${1}"cram_new.txt"
while read path
do
    cram_name=$(basename "$path" .cram)
    echo $cram_name

    dx download Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/${1}/"${cram_name}".cram -o cram/
    dx download Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/${1}/"${cram_name}".cram.crai -o cram/

    bash run.tryptase.sh $cram_name

    dx upload -r genotype/${cram_name} --path /Analysis_1/genotype/${cram_name}

    rm cram/${cram_name}.*
    #((i=i+1))
    #if [[ $i == 5 ]]
    #then
    #    break
    #fi
done < $input
