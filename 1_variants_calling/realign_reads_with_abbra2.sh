#!/bin/bash
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y # Error stream is merged with the standard output
#$ -l h_data=48G,h_rt=6:00:00
#$ -pe shared 1
#$ -r n # job is NOT rerunable
#$ -o log/realign.target.chrom.$TASK_ID.txt 
#$ -e log/realign.target.chrom.$TASK_ID.txt



set -euo pipefail

BAM=data_shortread/bam/MLDS00370PTAMPBAligned.sortedByCoord.out.bam
prefix=$(basename $BAM .bam)
OUT=data_shortread/bam/${prefix}.realigned_targets.2.bam
REF=GRCh38_no_alt_analysis_set_GCA_000001405.15.rmcontig.fasta
GTF=gencode.v32.annotation.gtf



date > log/realign.target.chrom.${SGE_TASK_ID}.txt 

java -Xmx16G -jar abra2-2.23.jar \
	    --in  ${BAM} \
	    --out ${OUT} \
	    --ref ${REF} \
	    --threads 4 \
	    --nosort \
	    --targets data_shortread/GATA1_alternative_regions.tmp.bed \
	    --gtf ${GTF} \
	    --dist 50000 --sua --tmpdir tmpdir --single 


echo "abra returns $?" 



