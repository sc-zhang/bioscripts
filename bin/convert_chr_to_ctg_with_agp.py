#!/usr/bin/env python
import sys


def reverse_seq(seq):
    base_db = {"A": "T", "T": "A", "G": "C", "C": "G"}
    rev_seq = ''.join([base_db[_] if _ in base_db else _ for _ in seq[::-1]])
    return rev_seq


def convert_chr_to_ctg(in_fa, in_agp, out_fa):
    print("Loading genome file")
    fa_db = {}
    with open(in_fa, 'r') as fin:
        for line in fin:
            if line[0] == '>':
                id = line.strip().split()[0][1:]
                fa_db[id] = []
            else:
                fa_db[id].append(line.strip().upper())

    for id in fa_db:
        fa_db[id] = ''.join(fa_db[id])

    print("Loading AGP file")
    ctg_db = {}
    with open(in_agp, 'r') as fin:
        for line in fin:
            if line.strip() == "" or line[0] == '#':
                continue
            data = line.strip().split()
            if data[4] != 'W':
                continue
            chrn = data[0]
            sp = int(data[1]) - 1
            ep = int(data[2])
            ctg = data[5]
            direct = data[-1]
            ctg_db[ctg] = [chrn, sp, ep, direct]

    print("Writing contig file")
    with open(out_fa, 'w') as fout:
        for ctg in sorted(ctg_db):
            chrn, sp, ep, direct = ctg_db[ctg]
            seq = fa_db[chrn][sp: ep]
            if direct == '-':
                seq = reverse_seq(seq)
            fout.write(">%s\n%s\n" % (ctg, seq))

    print("Finished")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python %s <in_fa> <in_agp> <out_fa>" % sys.argv[0])
    else:
        in_fa, in_agp, out_fa = sys.argv[1:]
        convert_chr_to_ctg(in_fa, in_agp, out_fa)
