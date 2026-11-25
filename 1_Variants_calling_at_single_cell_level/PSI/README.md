## Calculate PSI

Run the script `run_pysam_psi.umi_bc.copy.py` to calculate the PSI of regions of interest.
PSI = inclusion count / (inclusion count + exclusion count) normalized by length

Usage:
```
python run_pysam_psi.umi_bc.copy.py \
       -b file.bed \
       -m metadata.tsv \
       -o output.tsv \
```
Where the metadata file is a tab-separated file with the bam files you wish to use (Similar to the fofn file):

```    
bam   sample   cluster
sample1.cls1.bam   sample1       cluster1
sample1.cls2.bam   sample1       cluster2
```

The bed file should contain the regions of interest, with the following columns:

```
chrom	start	       end	       name	       score	strand
chrX	48783180	48783220	alt_tss       -1	+
chrX	48786585	48786590	tss           -1	+
chrX	48791852	48791980	intron        -1	+
```

The `name` column is used to identify the region in the output file.
The `score` column is not used, but it is required for compatibility. Fill in with -1 or 0.


