import pysam
import pandas as pd
from collections import defaultdict
import itertools
import HTSeq
import sys 
import pyranges as pr
import multiprocessing as mp
from glob import glob


def process_CIGAR( gPOS, cigar):
    block_curr = 0
    blocks = []

    for (ctype, length) in cigar:
        if ctype in (0, 7, 8):   # M,=,X        
            block_curr += length
        elif ctype == 2: #D
            block_curr += length
        elif ctype == 3: #N
            blocks.append((gPOS, gPOS + block_curr))
            gPOS += (block_curr + length)
            block_curr = 0

    blocks.append((gPOS, gPOS + block_curr))
    return blocks


def get_bc_info(metadata_file):

    bc = pd.read_csv(metadata_file, sep = ',')
    BC_CellType = defaultdict(lambda : defaultdict(lambda: 'NA'))

    for row in bc.itertuples():
        sample  = row._2
        barcode = row.cellbarcode.removeprefix(sample + "_")
        assert len(barcode) == 16 and all(c in "ACTGN" for c in barcode)
 
        celltype = row.cell_type
        
        BC_CellType[sample][barcode] = celltype

    return  BC_CellType


def extraxt_gtf_chains(gtf_file):

    tx_annot = defaultdict(list)
    gn_annot = defaultdict(lambda: defaultdict(lambda: 'NA'))

    gtf_file = pysam.TabixFile(gtf_file)
    for record in gtf_file.fetch(parser = pysam.asGTF()):
        if record.source == "HAVANA":
            
            attr_dict = {k: v for k, v in (x.strip().split(' ') for x in record.attributes.split(';') if x.strip())}
            gene_name = attr_dict['gene_name'].strip('"')

            if record.feature == "exon":
                exon_str = HTSeq.GenomicInterval(record.contig, record.start, record.end, record.strand)
                tx_annot[(record.transcript_id, record.strand)].append(exon_str)
            elif record.feature == "gene":
                gn_annot[gene_name]['coords'] = (record.contig, record.start, record.end)
            elif record.feature == "transcript":
                try: 
                    gn_annot[gene_name]['transcripts'].append((record.transcript_id, record.strand))
                except: 
                    gn_annot[gene_name]['transcripts'] = [(record.transcript_id, record.strand)]

    sys.stderr.write(f"Found {len(tx_annot)} canonical transcripts.\n")

    for (tx_id, tx_strand), tx_exons in tx_annot.items():

        tx_exons = sorted(tx_exons, key = lambda x: (x.start, x.end))

        introns_ss_5 = []
        introns_ss_3 = []

        for i in range(0, len(tx_exons) - 1):
            introns_ss_5.append(tx_exons[i].end)
            introns_ss_3.append(tx_exons[i + 1].start)

        tx_annot[(tx_id, tx_strand)] = [introns_ss_5, introns_ss_3]

    return tx_annot, gn_annot


if __name__ == "__main__":

    metadata_file = sys.argv[1]
    #BC_CellType = get_bc_info(metadata_file)
    BC_ReadID_dict = defaultdict(lambda : defaultdict(lambda: 'NA'))

    dataset, bam_file, gtf_file = sys.argv[2 : 5]

    tx_annot, gn_annot = extraxt_gtf_chains(gtf_file)

    sys.stderr.write("Finished loading transcript exons.\n")
    
    bh = pysam.AlignmentFile(bam_file , 'rb')

    for gene, gene_vals in gn_annot.items():

        sub_tx_annot = {k: v for k, v in tx_annot.items() if k in gene_vals['transcripts']}
        print(f"gene {gene} coords: {gene_vals['coords']}", flush=True)

        for j, read in enumerate(bh.fetch(*gene_vals['coords']), 1):

            if j % 5000 == 0:
                sys.stderr.write(f"\tProcessed {j:,} reads from dataset {dataset}.\n")

            read_name = read.query_name

            BC = read.get_tag('CB')

            cell_type = 'cell_type' #BC_CellType[dataset][BC]

            exons = process_CIGAR(read.positions[0], read.cigartuples)

            read_ss_5_tmp = []
            read_ss_3_tmp = []

            for jx in range(0, len(exons) - 1):
                read_ss_5_tmp.append(exons[jx][1])
                read_ss_3_tmp.append(exons[jx + 1][0])

            if read.is_reverse:
                read_ss_5 = read_ss_3_tmp[::-1]
                read_ss_3 = read_ss_5_tmp[::-1]
            else:
                read_ss_5 = read_ss_5_tmp
                read_ss_3 = read_ss_3_tmp

            read_chain_str = ':'.join([f"{ss5}-{ss3}" for ss5, ss3 in zip(read_ss_5, read_ss_3)])

            matches = [[101, False, 'NA']]

            if not read.is_reverse:
                continue


            for (tx_annot_id, tx_annot_strand), exons_annot in sub_tx_annot.items():

                if len(exons) >= 2:
                    if tx_annot_strand == "-":
                        exons_annot_ss_5 = exons_annot[1][::-1]
                        exons_annot_ss_3 = exons_annot[0][::-1]
                    else:
                        exons_annot_ss_5 = exons_annot[0]
                        exons_annot_ss_3 = exons_annot[1]

                    tx_chain_str = ':'.join([f"{ss5}-{ss3}" for ss5, ss3 in zip(exons_annot_ss_5, exons_annot_ss_3)])

                    delta_in = read_chain_str in tx_chain_str
                                    
                    delta = 0
                    for kx in range(0, len(read_ss_5)):

                            if kx >= len(exons_annot_ss_5):
                                delta += 1e7 ; break

                            delta += abs(exons_annot_ss_5[kx] - read_ss_5[kx])
                            delta += abs(exons_annot_ss_3[kx] - read_ss_3[kx])

                    matches.append([delta, delta_in, tx_annot_id])

            matches.sort(key = lambda x: x[0])
            print("log", read_name, read.is_reverse, dataset, len(exons), matches[0])

            ov_len, contained, tx_id_lab = matches[0]
            
            outline = ["out", dataset, read_name, BC, cell_type, tx_id_lab, ov_len, contained, gene]
            print('\t'.join(map(str, outline)))

    sys.stderr.write(f"Finished processing dataset {dataset}.\n")




        
