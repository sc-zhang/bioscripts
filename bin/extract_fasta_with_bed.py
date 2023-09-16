#!/usr/bin/env python
import sys


def rev_seq(seq):
    rseq = ""
    base_db = {"A": "T", "T": "A", "C": "G", "G": "C"}
    for base in seq[::-1]:
        if base in base_db:
            rseq += base_db[base]
        else:
            rseq += base
    return rseq


def extract_fa_with_bed(in_fa, in_bed, out_fa):
    print("Loading fasta")
    fa_db = {}
    with open(in_fa, 'r') as fin:
        seq = ""
        id = ""
        for line in fin:
            if line[0] == '>':
                if seq:
                    fa_db[id] = seq
                id = line.strip().split()[0][1:]
                seq = ""
            else:
                seq += line.strip().upper()
        if seq:
            fa_db[id] = seq
    
    print("Loading bed and writing fasta")
    with open(in_bed, 'r') as fin:
        with open(out_fa, 'w') as fout:
            for line in fin:
                data = line.strip().split()
                chrn = data[0]
                sp = int(data[1])-1
                ep = int(data[2])
                direct = data[3]
                id = data[4]
                seq = fa_db[chrn][sp: ep]
                if direct == '-':
                    seq = rev_seq(seq)
                fout.write(">%s\n%s\n"%(id, seq))
    print("Finsihed")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python %s <in_fa> <in_bed> <out_fa>"%sys.argv[0])
        print("Notice: bed should be 5 columns: \"ID, start, end, direction, id\", positions should be 1-based")
    else:
        in_fa, in_bed, out_fa = sys.argv[1:]
        extract_fa_with_bed(in_fa, in_bed, out_fa)
