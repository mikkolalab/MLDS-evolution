#!/bin/bash
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y # Error stream is merged with the standard output
#$ -l h_data=126G,h_rt=9:00:00 
#$ -pe shared 1
#$ -r n # job is NOT rerunable
#$ -o log/star.$TASK_ID.txt 
#$ -e log/star.$TASK_ID.txt


set -euo pipefail

# Load module
. /u/local/Modules/default/init/modules.sh
module load star


FO=data_shortread/fastq/MLDS00370PTAMPB_S169_L004_R1_001.fastq.gz.extracted.fastq.gz
FQ=data_shortread/fastq/MLDS00370PTAMPB_S169_L004_R2_001.fastq.gz.extracted.fastq.gz
WL=data_shortread/fastq/MLDS00370PTAMPB_S169_L004_R1_001.fastq.gz.whitelist.tsv
WLF=data_shortread/fastq/MLDS00370PTAMPB_S169_L004_R1_001.fastq.gz.whitelist.filtered.tsv

INDEX=/u/project/gxxiao/giovas/data/STAR_hg38_index
GTF=/u/project/gxxiao/giovas/data/hg38/gencode.v32.annotation.gtf


date > log/star.${SGE_TASK_ID}.txt



less ${WL} | cut -f1 > ${WLF}


STAR \
	--runThreadN 8 \
	--genomeDir ${INDEX} \
        --readFilesIn ${FQ} ${FO} \
	--readFilesCommand zcat \
	--soloType CB_UMI_Simple \
	--soloCBstart 1 \
	--soloCBlen 16 \
	--soloUMIstart 17 \
	--soloUMIlen 12 \
	--soloFeatures Gene \
	--soloCBwhitelist ${WLF} \
	--soloBarcodeReadLength 122 \
        --outReadsUnmapped None \
        --outSAMunmapped Within \
        --alignSJDBoverhangMin 6 \
        --alignIntronMax 500000 \
	--alignIntronMin 50 \
        --alignSplicedMateMapLminOverLmate 0 \
        --alignSplicedMateMapLmin 30 \
	--seedSearchStartLmax 12 \
	--alignSJoverhangMin 4 \
	--sjdbGTFfile ${GTF} \
	--twopassMode Basic \
        --outSAMattrRGline ID:GRPundef \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix data_shortread/bam/MLDS00370PTAMPB 



