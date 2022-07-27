#!/usr/bin/env python
import sys
import os


def split_fa(in_fa, in_agp, out_dir):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    fa_db = {}
    with open(in_fa, 'r') as fin:
        for line in fin:
            if line[0] == '>':
                id = line.strip().split()[0][1:]
                fa_db[id] = []
            else:
                fa_db[id].append(line.strip())

    for id in fa_db:
        fa_db[id] = ''.join(fa_db[id])

    chr_ctgs = {}
    with open(in_agp, 'r') as fin:
        for line in fin:
            if line.strip() == "" or line[0] == '#':
                continue
            data = line.strip().split()
            if data[4] != 'W':
                continue
            chrn = data[0]
            ctg = data[5]
            if chrn==ctg:
                chrn = 'Unanchored'
            if chrn not in chr_ctgs:
                chr_ctgs[chrn] = []
            chr_ctgs[chrn].append(ctg)

    for chrn in chr_ctgs:
        out_fn = os.path.join(out_dir, "%s.fasta"%chrn)
        with open(out_fn, 'w') as fout:
            for id in chr_ctgs[chrn]:
                fout.write(">%s\n%s\n"%(id, fa_db[id]))


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python %s <in_fa> <in_agp> <out_dir>"%sys.argv[0])
    else:
        in_fa, in_agp, out_dir = sys.argv[1:]
        split_fa(in_fa, in_agp, out_dir)
        
