#!/bin/bash

"********************************************************************* 
* Atlas Khan
* Kiryluk Lab (http://www.columbiamedicine.org/divisions/kiryluk/) 
* TPSABE1 project
* Version 1.0.0 
* (C) 2023 Nephrology Dep Medicine 
* Columbia University Medical Center
*********************************************************************"


#dir="/opt/notebooks/IGM_Joshua/Analysis_1"

dir="/opt/notebooks/Analysis_1"
echo "Notice current working directory!!!"
echo $dir

sample=$1

cfile="${dir}/cram/${sample}.cram"
odir="${dir}/genotype/${sample}"
mkdir -p ${odir}
tdir="${dir}/tmp/${sample}"
mkdir -p ${tdir} 
rdir="${dir}/reference" #### fasta file from Josh

bash ./tryptase.sh ${sample} ${cfile} ${odir} ${tdir} ${rdir}


