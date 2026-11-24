import pysam
import sys
from optparse import OptionParser
from collections import Counter



def getOptions():
    parser = OptionParser()
    parser.add_option("--f", dest = "infile", help = "Input file",
                      metavar = "FILE", type = "string", default = "")
    parser.add_option("--o", dest = "outfile",
                      help = "output file", metavar = "FILE", type = "string")
    (options, args) = parser.parse_args()
    return options



def reverseComplement(seq):
    """ Returns the reverse complement of a DNA sequence, 
        retaining the case of each letter"""
    complement = ""

    for base in seq:
        if base == "A": complement += "T"
        elif base == "T": complement += "A" 
        elif base == "G": complement += "C"
        elif base == "C": complement += "G"
        elif base == "N": complement += "N"
        elif base == "a": complement += "t"
        elif base == "t": complement += "a" 
        elif base == "g": complement += "c"
        elif base == "c": complement += "g"
        elif base == "n": complement += "n"
        elif base == "*": complement += "*"
        else:
            complement += base
            print("Warning: reverse complement function encountered unknown base " + "'" + base + "'")

    reverseComplement = complement[::-1]

    return reverseComplement

def main():

    polyA = "A"*20
    polyT = "T"*20

    options = getOptions() 

    bam_file_in  = pysam.AlignmentFile(options.infile, "rb", check_sq=False)
    bam_file_out = pysam.AlignmentFile(options.outfile, "wb", template = bam_file_in)

    i = 0
    counter = Counter()
    for record in bam_file_in:
        counter['all'] += 1

        if polyT in record.seq[0:50]:
            record.seq = reverseComplement(record.seq)
            counter['flipped'] += 1
            try: record.qual = record.qual[::-1]
            except: pass 

        if (not polyT in record.seq) and (not polyA in record.seq):
            counter['no poly A'] += 1
            continue


        bam_file_out.write(record)
        i += 1
        if i % 100000 == 0:
            sys.stderr.write(f"Processed {i:,} reads. Counter = {counter}\n")
            sys.stderr.flush()


if __name__ == '__main__':
    main()
