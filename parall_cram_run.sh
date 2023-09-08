#!/bin/bash

SECONDS=0
start=$(date +%s)
i=0
N=10
input="bulk"${1}"cram_new.txt"
while read path
do
        cram_name=$(basename "$path" .cram)
        echo $cram_name
    
        dx download Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/${1}/"${cram_name}".cram -o cram/
        dx download Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/${1}/"${cram_name}".cram.crai -o cram/
    (
        dir=$(pwd)

        sample=$cram_name

        cfile="${dir}/cram/${sample}.cram"
        odir="${dir}/genotype/${sample}"
        mkdir -p ${odir}
        tdir="${dir}/tmp/${sample}"
        mkdir -p ${tdir} 
        rdir="${dir}/reference" #### fasta file from Josh

        echo "running tryptase.sh"

        mkdir -p ${odir}
        mkdir -p ${tdir}

        samtools="/usr/local/bin/samtools"
        bwa="${dir}/bwa-master/bwa"
        python="${dir}/Python-2.7/python"
        region="chr16:1200000-1299999" # GRCh38/hg38
        fa_cons="${rdir}/CONS.fa"
        fa_tps="${rdir}/TPS.fa"

        bam_raw="${odir}/${sample}.tryptase.bam"
        fq_raw="${odir}/${sample}.tryptase.fq"
        bam_aln="${odir}/${sample}.cons.sam"
        fa_hap="${odir}/${sample}.fa"
        bam_hap="${odir}/${sample}.consx.sam"
        out_hap="${odir}/${sample}.output.txt"



        ## Extract reads mapping to the general tryptase locus from original BAM file
        ${samtools} view --no-PG -F 0xF0C -b -o ${bam_raw} ${cfile} ${region} 



        ## Convert to interleaved FASTQ
        ${samtools} collate -Ou ${bam_raw} ${tdir} | ${samtools} fastq > ${fq_raw}



        ## Map reads to consensus sequence, retaininly only mapped reads
        ${bwa} mem -p -M -R '@RG\tID:1\tSM:'${sample} ${fa_cons} ${fq_raw} | awk '$1 != "@PG"' | ${samtools} view --no-PG -F 0xF0C -S -h -o ${bam_aln}



        ## Cluster mapped reads into distinct haplotypes
        cat ${bam_aln} | awk '$10!~/N/' | ${python} parseHaplotypes.py > ${fa_hap}

        cat ${fa_hap} | \
        sed 's/ .*//' | \
        sed 's/^X*//' | \
        sed 's/X*$//' | \
        tr X N | \
        ${bwa} mem -M ${fa_tps} - | \
        awk '$1 != "@PG"' | \
        awk '$1!~/_[0-9]_/' > \
        ${bam_hap}

        cat ${bam_hap} | \
        grep -v ^@ | \
        cut -f1,3,4,6,12-13 | \
        tr _ "\t" | \
        awk '$3>5' | \
        sort -k4,4 -k3nr > \
        ${out_hap}

        dx upload -r genotype/${cram_name} --path /Analysis_1/genotype/${cram_name}

        $samtools depth cram/${cram_name}.cram | awk -v cram_var="$cram_name" 'OFS="\t" {sum+=$3} END { print cram_var,sum/NR}' > ${dir}/tsv/${cram_name}.sample.tsv

        awk -v cram_var="$cram_name" 'OFS="\t" {print cram_var,$0}' genotype/${cram_name}/${cram_name}.output.txt > ${dir}/tsv/${cram_name}.genotype.tsv

        rm cram/${cram_name}.*
    ) &

    # allow to execute up to $N jobs in parallel
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
        # now there are $N jobs already running, so wait here for any job
        # to be finished so there is a place to start next one.
        wait -n
    fi

    #after running certain number of jobs, terminate
    #((i=i+1))
    #if [[ $i == 100 ]]
    #then
    #   break
    #fi
done < $input

# no more jobs to be started but wait for pending jobs
# (all need to be finished)
wait

echo "all done"
duration=$SECONDS
echo "start: $start"
echo "end: $(date +%s)"
echo "duration: $(($duration / 3600)) hours, $((($duration / 60) % 60)) minutes and $(($duration % 60)) seconds elapsed."


