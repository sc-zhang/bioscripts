#!/usr/bin/env python
import argparse


def get_opt():
    group = argparse.ArgumentParser()
    group.add_argument('-i', '--input', help='Input maf file', required=True)
    group.add_argument('-o', '--output', help='Out put file', required=True)
    return group.parse_args()


def convert_anchorwave(in_file, out_file):
    print("Converting")
    with open(in_file, 'r') as fin:
        with open(out_file, 'w') as fout:
            tmp = []
            chrs = []
            fout.write("#Ref\tStart\tEnd\tQuery\tStart\tEnd\tType\n")
            for line in fin:
                if line.strip() == '' or line[0] == '#' or line[0] == 'a':
                    continue
                data = line.strip().split()
                chrs.append(data[1])
                tmp.append(data[6])
                if len(tmp) == 2:
                    seq_len = len(tmp[0])
                    print("\tConverting pair: %s, length: %d"%(','.join(chrs), seq_len))
                    ref_pos = 0
                    qry_pos = 0
                    var_info = []
                    per_cnt = int(seq_len / 10)
                    print("\t", end="")
                    for _ in range(seq_len):
                        if (_+1)%per_cnt == 0:
                            print("%d%%"%(int((_+1)/per_cnt)*10), end='\t', flush=True)
                        ref_base = tmp[0][_]
                        qry_base = tmp[1][_]
                        if ref_base == '-':
                            var_type = 'INS'
                            var_info.append([ref_pos, qry_pos, var_type])
                            qry_pos += 1
                        elif qry_base == '-':
                            var_type = 'DEL'
                            var_info.append([ref_pos, qry_pos, var_type])
                            ref_pos += 1
                        elif ref_base != qry_base:
                            var_type = 'SNP'
                            var_info.append([ref_pos, qry_pos, var_type])
                            ref_pos += 1
                            qry_pos += 1
                        else:
                            ref_pos += 1
                            qry_pos += 1
                    print()
                    if len(var_info) == 0:
                        tmp = []
                        chrs = []
                        continue
                    print("\tMerging pair: %s, length: %d"%(','.join(chrs), seq_len))
                    merge_info = [[var_info[0][0], var_info[0][0], 
                                var_info[0][1], var_info[0][1], 
                                var_info[0][2]]]
                    for _ in range(1, len(var_info)):
                        cur_info = var_info[_]
                        if cur_info[-1] == merge_info[-1][-1]:
                            is_continue = False
                            if cur_info[0] == merge_info[-1][1] + 1:
                                merge_info[-1][1] = cur_info[0]
                                is_continue = True
                            if cur_info[1] == merge_info[-1][3] + 1:
                                merge_info[-1][3] = cur_info[1]
                                is_continue = True
                            if not is_continue:
                                merge_info.append([cur_info[0], cur_info[0],
                                                cur_info[1], cur_info[1],
                                                cur_info[2]])
                        else:
                            merge_info.append([cur_info[0], cur_info[0],
                                            cur_info[1], cur_info[1],
                                            cur_info[2]])
                    print("\tWriting pair: %s, length: %d"%(','.join(chrs), seq_len))
                    for rsp, rep, qsp, qep, var_type in merge_info:
                        fout.write("%s\t%d\t%d\t%s\t%d\t%d\t%s\n"%(chrs[0], rsp+1, rep+1, 
                                                                   chrs[1], qsp+1, qep+1,
                                                                   var_type))
                    chrs = []
                    tmp = []
    print("Finished")


if __name__ == "__main__":
    opts = get_opt()
    in_file = opts.input
    out_file = opts.output
    convert_anchorwave(in_file, out_file)

