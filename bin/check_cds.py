#!/usr/bin/env python3
import argparse
from enum import Enum


class CDS_Type(Enum):
    VALID = 1
    LENGTH_ERROR = 2
    MISSING_START_CODON = 3
    MISSING_STOP_CODON = 4
    EARLY_STOP_CODON = 5


def get_opts():
    group = argparse.ArgumentParser()
    group.add_argument("-i", "--input", help="Input CDS file", required=True)
    group.add_argument(
        "--detail", help="If set, output detail information", action="store_true"
    )
    group.add_argument(
        "-o",
        "--output",
        help="Output summary file, if not set, output to stdout",
        default="",
    )
    return group.parse_args()


def check_cds(in_cds, is_detail, out_summary):
    start_codon = set(["ATG"])
    stop_codon = set(["TAG", "TAA", "TGA"])
    cds_db = {}
    with open(in_cds, "r") as fin:
        for line in fin:
            if line.strip() == "":
                continue
            if line[0] == ">":
                gid = line.strip().split()[0][1:]
                cds_db[gid] = []
            else:
                cds_db[gid].append(line.strip().upper())

    for gid in cds_db:
        cds_db[gid] = "".join(cds_db[gid])

    detail_db = {}
    for gid in cds_db:
        detail_db[gid] = CDS_Type.VALID
        if len(cds_db[gid]) % 3 != 0:
            detail_db[gid] = CDS_Type.LENGTH_ERROR
        else:
            if cds_db[gid][:3] not in start_codon:
                detail_db[gid] = CDS_Type.MISSING_START_CODON
            elif cds_db[gid][-3:] not in stop_codon:
                detail_db[gid] = CDS_Type.MISSING_STOP_CODON
            else:
                for _ in range(3, len(cds_db[gid]) - 3, 3):
                    if cds_db[gid][_ : _ + 3] in stop_codon:
                        detail_db[gid] = CDS_Type.EARLY_STOP_CODON

    # Valid, length error, missing start codon, missing stop codon, early stop codon
    summary_info = [0, 0, 0, 0, 0]
    for gid in detail_db:
        summary_info[detail_db[gid].value - 1] += 1

    out_info = []
    out_info.append("# Summary")
    out_info.append("Valid:               %d" % summary_info[0])
    out_info.append("Length error:        %d" % summary_info[1])
    out_info.append("Missing start codon: %d" % summary_info[2])
    out_info.append("Missing stop codon:  %d" % summary_info[3])
    out_info.append("Early stop codon:    %d" % summary_info[4])

    if is_detail:
        out_info.append("")
        out_info.append("# Error detail")
        for gid in sorted(detail_db):
            if detail_db[gid] != CDS_Type.VALID:
                out_info.append("%s: %s" % (gid, detail_db[gid].name))

    if out_summary:
        with open(out_summary, "w") as fout:
            fout.write("%s\n" % ("\n".join(out_info)))
    else:
        print("%s" % ("\n".join(out_info)))


if __name__ == "__main__":
    opts = get_opts()
    in_cds = opts.input
    is_detail = True if opts.detail else False
    out_summary = opts.output
    check_cds(in_cds, is_detail, out_summary)
