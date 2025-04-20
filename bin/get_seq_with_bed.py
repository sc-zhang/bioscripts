#!/usr/bin/env python
import sys
import gzip


def reverse_seq(seq):
    base_db = {"A": "T", "T": "A", "C": "G", "G": "C"}
    rev_seq = "".join([base_db[_] if _ in base_db else _ for _ in seq.upper()[::-1]])
    return rev_seq


def get_seq_with_bed(in_fa, in_bed, out_fa):
    print("Loading bed")
    bed_db = {}
    with open(in_bed, "r") as fin:
        for line in fin:
            data = line.strip().split()
            chrn = data[0]
            sp = int(data[1]) - 1
            ep = int(data[2])
            if len(data) > 4:
                direct = data[3]
            else:
                direct = "+"
            gid = data[-1]
            if chrn not in bed_db:
                bed_db[chrn] = []
            bed_db[chrn].append([sp, ep, direct, gid])

    print("Extracting")
    if in_fa.endswith(".gz"):
        fin = gzip.open(in_fa, "rt")
    else:
        fin = open(in_fa, "r")

    fa_db = {}
    for line in fin:
        if line[0] == ">":
            cid = line.strip()[1:]
            fa_db[cid] = []
        else:
            fa_db[cid].append(line.strip())

    fin.close()
    for _ in fa_db:
        fa_db[_] = "".join(fa_db[_])

    with open(out_fa, "w") as fout:
        for cid in bed_db:
            for sp, ep, direct, gid in bed_db[cid]:
                fout.write(
                    ">%s\n%s\n"
                    % (
                        gid,
                        (
                            fa_db[cid][sp:ep]
                            if direct == "+"
                            else reverse_seq(fa_db[cid][sp:ep])
                        ),
                    )
                )

    print("Finished")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python %s <in_fa> <in_bed> <out_fa>" % sys.argv[0])
    else:
        in_fa, in_bed, out_fa = sys.argv[1:]
        get_seq_with_bed(in_fa, in_bed, out_fa)
