#!/bin/bash
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y # Error stream is merged with the standard output
#$ -l h_data=58G,h_rt=6:00:00
#$ -pe shared 1 
#$ -r n # job is NOT rerunable
#$ -o log/bam_bc_ds.$TASK_ID.txt 
#$ -e log/bam_bc_ds.$TASK_ID.txt


ja=job.lr_bam.txt 


PARMS=($(awk "NR==$SGE_TASK_ID" $ja))
bam=${PARMS[0]}


date > log/bam_bc_ds.${SGE_TASK_ID}.txt

python TALON/barcode_downsample_bam.py ${bam}

echo "returns $?" 

