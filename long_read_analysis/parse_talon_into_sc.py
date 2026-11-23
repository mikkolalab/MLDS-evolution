import os
import sys
import logging
from scipy import sparse
from scipy.io import mmwrite
from glob import glob
import pandas as pd
import numpy as np
from collections import Counter, defaultdict
import argparse





if __name__ == "__main__":

	logging.basicConfig(level = logging.INFO, format =  "%(asctime)s - %(levelname)s - %(message)s", force = True)
	logging.info("Starting processing...")

	parser = argparse.ArgumentParser(description = "Parse TALON output into single-cell format")
	parser.add_argument("sample", type = str, help = "Sample name to process")
	parser.add_argument("out_dir", type = str, help = "Output file out_dir (without extension)")
	args = parser.parse_args()

	sample = args.sample

	if not os.path.exists(args.out_dir):
		os.makedirs(args.out_dir)

	out_file_matrix   = os.path.join(args.out_dir, "matrix.mtx")
	out_file_barcodes = os.path.join(args.out_dir, "barcodes.tsv")
	out_file_genes    = os.path.join(args.out_dir, "genes.tsv")


	files = glob(f'data_longread/TALON_OUT/{sample}.tagged.p*_out_talon_read_annot.tsv')

	logging.info(f"Found {len(files)} files to process.")


	counter = defaultdict(dict)
	gata_counter = Counter()

	for read_annot in files:
		read_annot_df = pd.read_csv(read_annot, sep = '\t')
		for row in read_annot_df.itertuples():

			if row.annot_transcript_id.startswith("novel"):
				continue

			try: 
				counter[row.dataset][row.annot_transcript_id] += 1 
			except: 
				counter[row.dataset][row.annot_transcript_id] = 1

			if row.annot_gene_name == "GATA1":
				gata_counter[row.annot_transcript_id] += 1

		logging.info(f"Processed {read_annot}")



	for gata_tx, count in gata_counter.items():
		logging.info(f"GATA1_transcript: {sample} {gata_tx} {count}")

	logging.info("Finished processing files. Now writing output...")

	mat = pd.DataFrame(counter).fillna(0).astype(int)
	mat = mat.sort_index(axis = 1)

	counts = sparse.csr_matrix(mat.values) 
	mmwrite(out_file_matrix, counts)

	# write column names in order they appear
	with open(out_file_barcodes, 'w') as f:
		for col in mat.columns:
			f.write(col + '\n')

	# write row names in order they appear
	with open(out_file_genes, 'w') as f:
		for row in mat.index:
			f.write(f'{row}\t{row}\n')
