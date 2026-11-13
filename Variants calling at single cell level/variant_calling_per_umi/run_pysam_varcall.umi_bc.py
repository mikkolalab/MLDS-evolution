import pandas as pd
import pysam
from collections import defaultdict, Counter
import os
import logging
import argparse

rev = {"A": "T", "T": "A", "C": "G", "G": "C"}

if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(description="Run variant calling on BAM files using pysam")
    parser.add_argument('-v', '--variants', type=str, required=True, help='Path to the variants file (TSV format)')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output file name for the variants summary')
    parser.add_argument('-m', '--metadata', type=str, required=True, help='Path to the metadata file (TSV format)')
    args = parser.parse_args()

    VARS = pd.read_csv(args.variants, sep=",")
    VARS.columns = VARS.columns.str.strip().str.lower()

    metadata = pd.read_csv(args.metadata, sep="\t")
    metadata.columns = metadata.columns.str.strip().str.lower()

    logging.info(f"Starting variant calling: {VARS.shape[0]} variants")

    var_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    for i, row in metadata.iterrows():

        bam_file, sample, data_type = row['bam'], row['sm'], row['data_type']

        if not os.path.exists(bam_file + ".bai"):
            logging.error(f"Index file not found for {bam_file}")
            continue

        try:
            bh = pysam.AlignmentFile(bam_file)
        except Exception as e:
            logging.error(f"Error opening {bam_file}: {e}")
            continue

        for var in VARS.itertuples(index=False):
            var_annot = f"{var.chrom}_{var.pos}_{var.ref}_{var.alt}_{var.mut}"

            pileup = bh.pileup(var.chrom, var.pos - 2, var.pos + 1, 
                                stepper = "nofilter", truncate = True, max_depth = int(1e8))

            for pileupcolumn in pileup:
                if pileupcolumn.pos == (var.pos - 1):
                    for pileupread in pileupcolumn.pileups:
    
                        read_id = pileupread.alignment.query_name

                        if pileupread.alignment.is_secondary:                        
                            continue 
                        elif pileupread.is_refskip:
                            continue
 
                        try:
                            umi = pileupread.alignment.get_tag('UB')
                            bc  = pileupread.alignment.get_tag('CB')
                        except:
                            parts = read_id.split("_")
                            if len(parts) >= 3:
                                read_name, bc, umi = parts[0], parts[1], parts[2]
                            else:
                                continue

                        if len(var.ref) == 1 and len(var.alt) == 1:
                            if pileupread.is_del or (pileupread.query_position is None):
                                continue
                            base = pileupread.alignment.query_sequence[pileupread.query_position]

                        elif len(var.ref) > len(var.alt):
                            # deletion
                            if pileupread.indel == len(var.alt) - len(var.ref):
                                base = var.alt
                            else:
                                base = var.ref

                        elif len(var.ref) < len(var.alt):
                            # insertion
                            if pileupread.query_position is None:
                                continue
                            if pileupread.indel == len(var.alt) - len(var.ref):
                                base = var.alt
                            else:
                                base = var.ref
                        else:
                            base = var.ref  # 兜底

                        var_dict[(sample, bc, var_annot)][data_type][umi].append(base)

        bh.close()
        logging.info(f"\tDone processing {sample} ({i + 1}/{metadata.shape[0]})")

    outline = defaultdict(dict)

    for (sm, bc, var), val in var_dict.items():
        CHROM, POS, REF, ALT = var.split("_")[:4]

        for data_type, umi_bases in val.items():
            base_counter = Counter()
            for umi, bases in umi_bases.items():
                umi_base, count = Counter(bases).most_common(1)[0]
                if ALT in bases:
                    base_counter[ALT] += 1
                else:
                    base_counter[umi_base] += 1

            DP = sum(base_counter.values())
            AB = base_counter[ALT] / DP if DP > 0 else 0

            outline[(sm, bc)][f"{var}_{data_type}_cov"] = DP
            outline[(sm, bc)][f"{var}_{data_type}_ratio"] = AB

    table = pd.DataFrame(outline).T
    table = table.fillna('NA')
    table.to_csv(args.output, sep="\t")

    print(table)
