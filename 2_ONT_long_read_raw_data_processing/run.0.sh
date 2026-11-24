
REF='/u/project/gxxiao/giovas/data/hg38/GCF_000001405.40_GRCh38.p14_genomic.canonical_chroms.fa'
GTF='/u/project/gxxiao/giovas/data/hg38/gencode.v32.annotation.gtf'
TCDir='/u/project/gxxiao/giovas/projects/MLLT3_longread_project/scripts/TranscriptClean-2.1'


JXN=data/gencode.v32.junctions.tsv


rm -rf data/*


python ${TCDir}/accessory_scripts/get_SJs_from_gtf.py \
                        --f ${GTF} \
                        --g ${REF} \
                        --o ${JXN}


