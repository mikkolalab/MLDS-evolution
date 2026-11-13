#!/bin/bash
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y # Error stream is merged with the standard output
#$ -l h_data=24G,h_rt=1:00:00
#$ -pe shared 1
#$ -r n # job is NOT rerunable
#$ -o log/manual_quant.$TASK_ID.txt 
#$ -e log/manual_quant.$TASK_ID.txt



set -euo pipefail


#ja=job.bam_long_read.txt
ja=job.mecom.txt

PARMS=($(awk "NR==$SGE_TASK_ID" $ja))
SM=${PARMS[0]}
BAM=${PARMS[1]}

#GTF=gencode.v44.annotation.mllt3-s2.MLLT3_only.gtf.gz
GTF=gencode.v44.annotation.MLLT3_and_MECOM.gtf.gz 
METADATA_FILE=data_longread/orig_ident_cellbarcode.csv

date > log/manual_quant.${SGE_TASK_ID}.txt

# samtools index -@ 16 ${BAM}
# samtools view  -@ 16 -bh ${BAM} chr9:20341660-20622500 > ${BAM%.bam}_MLLT3.bam
# samtools index -@ 16 ${BAM%.bam}_MLLT3.bam
# echo "samtools returns $?"


python scripts/manual_read_assignment_to_isoforms.py \
    ${METADATA_FILE} \
    ${SM} \
    ${BAM} \
    ${GTF}

echo "python returns $?"
