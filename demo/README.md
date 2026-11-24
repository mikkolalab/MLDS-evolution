# Single-cell data processing



FQO=MLDS00272Pd0_S13_L006_R1_001.fastq.downsapled.gz
FQT=MLDS00272Pd0_S13_L006_R2_001.fastq.downsapled.gz



##

`
umi_tools whitelist \
      --stdin ${FQO} \
      --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNNN \
      --plot-prefix=${FQO}.whitelist \
      --subset-reads=200000000 \
      --log2stderr > ${FQO}.whitelist.tsv  
`

##

`
umi_tools extract \
    --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNNN \
    -I ${FQO} \
    -S ${FQO}.extracted.fastq.gz  \
    --read2-in=${FQT} \
    --read2-out=${FQT}.extracted.fastq.gz \
    --whitelist=${FQO}.whitelist.tsv 
`

## Align with STAR



FO=data_shortread/fastq/MLDS00370PTAMPB_S169_L004_R1_001.fastq.gz.extracted.fastq.gz
FQ=data_shortread/fastq/MLDS00370PTAMPB_S169_L004_R2_001.fastq.gz.extracted.fastq.gz
WL=data_shortread/fastq/MLDS00370PTAMPB_S169_L004_R1_001.fastq.gz.whitelist.tsv
WLF=data_shortread/fastq/MLDS00370PTAMPB_S169_L004_R1_001.fastq.gz.whitelist.filtered.tsv

INDEX=/u/project/gxxiao/giovas/data/STAR_hg38_index
GTF=/u/project/gxxiao/giovas/data/hg38/gencode.v32.annotation.gtf


date > log/star.${SGE_TASK_ID}.txt



less ${WL} | cut -f1 > ${WLF}

`
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

`

