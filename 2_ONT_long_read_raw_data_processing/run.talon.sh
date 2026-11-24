#!/bin/bash
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y # Error stream is merged with the standard output
#$ -l h_data=63G,h_rt=11:55:00
#$ -pe shared 1
#$ -r n # job is NOT rerunable
#$ -o log/talon_talon.$TASK_ID.txt 
#$ -e log/talon_talon.$TASK_ID.txt



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

GTF='/u/project/gxxiao/giovas/data/hg38/gencode.v32.annotation.gtf'



set -euo pipefail


date > log/talon_talon.${SGE_TASK_ID}.txt

echo $SM 

config=data_longread/TALON_OUT/${SM}.config
printf "${SM},ONT,data_longread/TALON_OUT/${SM}_lab_labeled.sam\n" > ${config}
cat ${config}


bam=data_longread/TALON_OUT/${SM}_lab_labeled.bam
if [ -f "$bam" ]
then
	samtools view -h $bam > data_longread/TALON_OUT/${SM}_lab_labeled.sam
fi



rm -f  data_longread/TALON_OUT/${SM}_talon.db


talon_initialize_database \
	--f ${GTF} \
        --g 'GRCh38' \
        --a 'v32' \
        --5p 500 \
        --3p 300 \
        --idprefix 'novel' \
        --o data_longread/TALON_OUT/${SM}_talon &> log/talon_talon.${SGE_TASK_ID}.err 
echo "LOG talon init returns $?"

talon \
        --f ${config} \
        --db data_longread/TALON_OUT/${SM}_talon.db \
	--cb \
        --build 'GRCh38' \
        --threads 12\
        --cov 0.5 \
	--identity 0.7 \
        --tmpDir data_longread/TALON_OUT/${SM}_talontmpdir \
        --o data_longread/TALON_OUT/${SM}_out

echo "LOG talon run returns $?"
rm -rf data_longread/TALON_OUT/${SM}_talontmpdir


talon_fetch_reads \
	--db data_longread/TALON_OUT/${SM}_talon.db \
	--build 'GRCh38' \
	--o data_longread/TALON_OUT/${SM}_ann
echo "LOG talon fetch returns $?"


samtools view -bSh -@ 8 data_longread/TALON_OUT/${SM}_lab_labeled.sam > data_longread/TALON_OUT/${SM}_lab_labeled.bam
echo "LOG bam conversion returns $?"


du -sh data_longread/TALON_OUT/${SM}_lab_labeled.?am
#rm -f  data_longread/TALON_OUT/${SM}_lab_labeled.sam


