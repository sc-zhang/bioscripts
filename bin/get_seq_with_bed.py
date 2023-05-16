#!/usr/bin/env python
import sys


def reverse_seq(seq):
    base_db = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    rev_seq = ''.join([base_db[_] if _ in base_db else _ for _ in seq.upper()[::-1]])
    return rev_seq


def get_seq_with_bed(in_fa, in_bed, out_fa):
    print("Loading bed")
    bed_db = {}
    with open(in_bed, 'r') as fin:
        for line in fin:
            data = line.strip().split()
            chrn = data[0]
            sp = int(data[1])-1
            ep = int(data[2])
            if len(data) > 4:
                direct = data[3]
            else:
                direct = '+'
            gid = data[-1]
            if chrn not in bed_db:
                bed_db[chrn] = []
            bed_db[chrn].append([sp, ep, direct, gid])
    
    print("Extracting")
    with open(in_fa, 'r') as fin:
        with open(out_fa, 'w') as fout:
            cid = ""
            seq = ""
            with open(in_fa, 'r') as fin:
                for line in fin:
                    if line.strip() == "":
                        continue
                    if line[0] == '>':
                        if seq != "":
                            for sp, ep, direct, gid in bed_db[cid]:
                                fout.write(">%s\n%s\n"%(gid, 
                                                        seq[sp: ep] if direct=='+' 
                                                        else reverse_seq(seq[sp: ep])))
                        cid = line.strip()[1:]
                        seq = ""
                    else:
                        seq += line.strip()
            if seq != "":
                for sp, ep, direct, gid in bed_db[cid]:
                                fout.write(">%s\n%s\n"%(gid, 
                                                        seq[sp: ep] if direct=='+' 
                                                        else reverse_seq(seq[sp: ep])))
    print("Finished")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python %s <in_fa> <in_bed> <out_fa>"%sys.argv[0])
    else:
        in_fa, in_bed, out_fa = sys.argv[1:]
        get_seq_with_bed(in_fa, in_bed, out_fa)
