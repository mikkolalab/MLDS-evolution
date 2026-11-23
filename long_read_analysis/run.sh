#!/bin/bash
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y # Error stream is merged with the standard output
#$ -l h_data=36G,h_rt=6:45:00
#$ -pe shared 1
#$ -r n # job is NOT rerunable
#$ -o log/talon_bc.$TASK_ID.txt 
#$ -e log/talon_bc.$TASK_ID.txt


# Load module
. /u/local/Modules/default/init/modules.sh
module load python/3.7.3


# Activate conda env
source activate base
conda activate isoquant



ja=job.lr_bam_bc.txt 
PARMS=($(awk "NR==$SGE_TASK_ID" $ja))
BAM=${PARMS[0]}

SM=$(basename $BAM ".bam")
REF='/u/project/gxxiao/giovas/data/hg38/GCF_000001405.40_GRCh38.p14_genomic.canonical_chroms.fa'
GTF='/u/project/gxxiao/giovas/data/hg38/gencode.v32.annotation.gtf'
TCDir='/u/project/gxxiao/giovas/projects/MLLT3_longread_project/scripts/TranscriptClean-2.1'
JXN='/u/project/gxxiao/giovas/projects/MLDS_project/TALON/data/gencode.v32.junctions.tsv'
DB='/u/project/gxxiao/giovas/projects/MLDS_project/TALON/data/gencode.v32.db'



set -euo pipefail


date > log/talon_bc.${SGE_TASK_ID}.txt
echo $SM 


python TALON/flip_reads.py \
	--f ${BAM} \
        --o ${BAM%'.bam'}.flipped.bam
echo "LOG flip reads returns $?" 


chroms=$(less ${REF}.fai | cut -f1)

samtools index -@ 4 ${BAM%'.bam'}.flipped.bam 
samtools view -h ${BAM%'.bam'}.flipped.bam $chroms > ${BAM%'.bam'}.refined.sam 
echo "LOG samtools returns $?"

rm -f ${BAM%'.bam'}.flipped.bam ${BAM%'.bam'}.flipped.bam.bai 



mkdir -p data_longread/TALON_OUT/${SM}_tctmp

python ${TCDir}/TranscriptClean.py \
        --correctMismatches False \
        --correctIndels False \
        --threads 24 \
        --sam ${BAM%'.bam'}.refined.sam \
        --genome ${REF} \
        --spliceJns ${JXN} \
        --outprefix data_longread/TALON_OUT/${SM}_tc \
        --maxSJOffset 5 \
        --tmpDir data_longread/TALON_OUT/${SM}_tctmp


echo "LOG transcriptclean returns $?"
rm -rf data_longread/TALON_OUT/${SM}_tctmp

mkdir -p data_longread/TALON_OUT/${SM}_labtmp

talon_label_reads \
        --f data_longread/TALON_OUT/${SM}*_tc_clean.sam \
        --g ${REF} \
        --ar 20 \
        --o data_longread/TALON_OUT/${SM}_lab \
        --tmpDir data_longread/TALON_OUT/${SM}_labtmp

echo "LOG talon labeling returns $?"
rm -rf data_longread/TALON_OUT/${SM}_labtmp

rm -f ${BAM%'.bam'}.refined.sam
rm -f data_longread/TALON_OUT/${SM}*_tc_clean.sam
rm -f data_longread/TALON_OUT/${SM}*.fa

