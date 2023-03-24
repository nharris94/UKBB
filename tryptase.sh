# ascertain tryptase content from sequence data

sample=${1}
cfile=${2}
odir=${3}
tdir=${4}
rdir=${5}

mkdir -p ${odir}
mkdir -p ${tdir}

samtools="/usr/local/bin/samtools"
bwa="/opt/notebooks/Analysis_1/bwa/bwa"
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
cat ${bam_aln} | awk '$10!~/N/' | /opt/notebooks/Analysis_1/Python-2.7/python parseHaplotypes.py > ${fa_hap}

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


