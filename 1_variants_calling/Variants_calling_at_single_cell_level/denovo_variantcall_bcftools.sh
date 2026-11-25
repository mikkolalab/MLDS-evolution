#!/bin/bash
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y # Error stream is merged with the standard output
#$ -l h_data=16G,h_rt=2:00:00,highp 
#$ -pe shared 1
#$ -r n # job is NOT rerunable
#$ -o SGE_STDOUT.txt 


ja=bam.job
PARMS=($(awk "NR==$SGE_TASK_ID" $ja))
BAM=${PARMS[0]}


. /u/local/Modules/default/init/modules.sh
module load bcftools
module load samtools




REF=${HOME}/gxxiao/data/GCF_000001405.40_GRCh38.p14_genomic.fa 

log=log/bcftools.${SGE_TASK_ID}.txt


declare -A genes 

genes=( ["chrX:48,783,671-48,801,157"]="GATA1"  )

genes=( ["chrX:48783671-48801157"]="GATA1"  )

date > $log 

PREFIX=${BAM%".bam"}.tmp

for coords in "${!genes[@]}"; 
do 

    GENE=${genes[$coords]}

    tmp=${PREFIX}.${GENE}.tmp.bam 

    samtools view -bh ${BAM} ${coords} > ${tmp}

    bcftools mpileup ${tmp} -f ${REF} | bcftools call -mv -Oz > ${PREFIX}.${GENE}.vcf.gz 2>> $log 

    bcftools index ${PREFIX}.${GENE}.vcf.gz 2>> $log

    rm -f $tmp 

    echo "returns $GENE $?" >> $log

done


 
