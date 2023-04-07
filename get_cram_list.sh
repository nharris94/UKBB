#!/bin/bash

find /mnt/project/Bulk/'Exome sequences'/'Exome OQFE CRAM files'/${1} -name *.cram > bulk10cram.txt

grep -vf <(ls /mnt/project/Analysis_1/genotype/) bulk${1}cram.txt > bulk${1}cram_new.txt