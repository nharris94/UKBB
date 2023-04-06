#!/bin/bash

i=0
N=10
input="bulk"${1}"cram_new.txt"
while read path
do
    (
        cram_name=$(basename "$path" .cram)
        echo $cram_name

        dx download Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/${1}/"${cram_name}".cram -o cram/
        dx download Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/${1}/"${cram_name}".cram.crai -o cram/

        bash run.tryptase.sh $cram_name

        # dx upload -r genotype/${cram_name} --path /Analysis_1/genotype/${cram_name}

        # rm cram/${cram_name}.*
    ) &

    # allow to execute up to $N jobs in parallel
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
        # now there are $N jobs already running, so wait here for any job
        # to be finished so there is a place to start next one.
        wait -n
    fi

    #after running certain number of jobs, terminate
    ((i=i+1))
    if [[ $i == 20 ]]
    then
       break
    fi
done < $input