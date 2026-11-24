#!/bin/bash
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y # Error stream is merged with the standard output
#$ -l h_data=36G,h_rt=4:00:00
#$ -pe shared 1
#$ -r n # job is NOT rerunable
#$ -o log/talon_matrix.$TASK_ID.txt 
#$ -e log/talon_matrix.$TASK_ID.txt
#$ -t 34-40 



ja=job.lr_bam.txt
PARMS=($(awk "NR==$SGE_TASK_ID" $ja))
BAM=${PARMS[0]}

SM=$(basename $BAM ".tagged.bam")
OUT=data_longread/TALON_OUT/${SM}_gene_bc_matrices


set -euo pipefail


date > log/talon_matrix.${SGE_TASK_ID}.txt

echo $SM 

mkdir -p ${OUT}

python TALON/parse_talon_into_sc.py ${SM} ${OUT}

echo  "parser returns $?"

