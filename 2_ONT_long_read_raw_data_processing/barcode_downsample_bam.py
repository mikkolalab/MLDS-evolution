import pysam
import sys
from collections import defaultdict
import multiprocess as mp


BC_SIZE = 1000

def split_bam(j, barcode_set):

    bh = pysam.AlignmentFile(input_bam, 'rb')

    bamf = input_bam.replace(".bam", f".p{j}.bam")
    bamh = pysam.AlignmentFile(bamf, 'wb', template = bh)


    for i, read in enumerate(bh.fetch(until_eof = True), 1):
        barcode = read.get_tag('CB')

        if barcode in barcode_set:
            bamh.write(read)

        if i % 500000 == 0:
            sys.stderr.write(f'\t{j} {i:,} records processed.\n')

    bh.close()
    bamh.close()
    sys.stderr.write(f'\tDone writting bam to {bamf}\n')





if __name__ == "__main__":

    global input_bam

    input_bam = sys.argv[1]
    
    bh = pysam.AlignmentFile(input_bam, 'rb')

    barcode_lst = set()

    for i, read in enumerate(bh.fetch(until_eof = True), 1):
        barcode = read.get_tag('CB')
        barcode_lst.add(barcode)
        if i % 500000 == 0:
            sys.stderr.write(f'\t{i:,} records processed.\n')

    bh.close()

    sys.stderr.write(f'{len(barcode_lst)} barcodes found.\n')
    barcode_lst = list(barcode_lst)


    # split barcodes into chunks
    barcode_chunks = [ set(barcode_lst[i : i + BC_SIZE]) for i in range(0, len(barcode_lst), BC_SIZE) ]
    
    n = len(barcode_chunks)

    sys.stderr.write(f'{n} workers\n')

    pool = mp.Pool(processes = n)

    jobs = []
    for j, bcset in enumerate(barcode_chunks):
        jobs.append(pool.apply_async(split_bam, (j, bcset)))

    pool.close()
    pool.join()

    sys.stderr.write("All jobs finished.\n")

    

       
