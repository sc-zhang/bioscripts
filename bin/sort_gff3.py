#!/usr/bin/env python
import sys
import re


def generate_new_id_by_gff(chr_pre, in_gff):
    id_db = {}
    gff_db = {}
    is_first_id = True
    with open(in_gff, 'r') as f_in:
        for line in f_in:
            data = line.strip().split()
            if line[0] == '#' or len(data) < 9:
                continue
            if data[2] == 'gene':
                chrn = data[0]
                s_p = int(data[3])
                id = data[8].split(";")[0].split("=")[1]
                if is_first_id:
                    print("Check ID: %s" % id)
                    is_first_id = False
                if chrn not in id_db:
                    id_db[chrn] = []
                id_db[chrn].append([s_p, id])
                gff_db[id] = []
            gff_db[id].append(line)

    new_id_db = {}
    ordered_id = []
    tig_base = 10
    for chrn in sorted(id_db, key=lambda x: int(re.findall(r'\d+', x)[0]) if len(re.findall(r'\d+', x)) > 0 else 1000):
        base = 10
        for info in sorted(id_db[chrn]):
            if chrn[:3].lower() == 'chr':
                idx, hap = re.findall(r"(\d+)([A-Z]*)", chrn)[0]
                idx = int(idx)
                if not hap:
                    hap = "G"
                new_id = chr_pre + ".%02d%s%07d" % (idx, hap, base)
                base += 10
            else:
                new_id = chr_pre + ".%08d" % tig_base
                tig_base += 10
            ordered_id.append(info[1])
            new_id_db[info[1]] = new_id
    return new_id_db, gff_db, ordered_id


def rename_id(chr_pre, in_gff, out_gff):
    print("Generating rename list")
    rename_id_db, gff_db, ordered_id = generate_new_id_by_gff(chr_pre, in_gff)

    print("Dealing gff3")
    with open(out_gff, 'w') as fout:
        fout.write("###gff version 3\n")
        for ori_id in ordered_id:
            mrna_idx = 1
            for line in gff_db[ori_id]:
                data = line.strip().split()
                if data[2] == 'gene':
                    gid = rename_id_db[ori_id]
                    data[8] = "ID=%s;Name=%s" % (gid, gid)
                elif data[2] == 'mRNA':
                    mrid = "%s.t%d" % (gid, mrna_idx)
                    mrna_idx += 1
                    other_idx_db = {}
                    data[8] = "ID=%s;Name=%s;Parent=%s" % (mrid, mrid, gid)
                else:
                    feature = data[2]
                    if feature not in other_idx_db:
                        other_idx_db[feature] = 1
                    other_id = "%s.%s%d" % (mrid, feature, other_idx_db[feature])
                    other_idx_db[feature] += 1
                    data[8] = "ID=%s;Name=%s;Parent=%s" % (other_id, other_id, mrid)
                fout.write("%s\n" % ("\t".join(data)))
            fout.write("###\n")
    print("Finished")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python " + sys.argv[0] + " <chr_prefix> <in_gff3> <out_gff3>")
        print("Notice: sort and rename id with in_gff by coordinate, the chromosome ID should be like: Chr01 for mono "
              "assembly, Chr01A for phased assembly.")
        print("Example: python " + sys.argv[0] + " CB5 in.gff3 out.gff3")
    else:
        chr_pre, in_gff3, out_gff3 = sys.argv[1:]
        rename_id(chr_pre, in_gff3, out_gff3)
