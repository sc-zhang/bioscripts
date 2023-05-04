#!/usr/bin/env python3
import argparse
import numpy as np
import bisect


def get_opts():
    group = argparse.ArgumentParser()
    group.add_argument('-r', '--ref', help='ref.bed', required=True)
    group.add_argument('-q', '--qry' ,help='qry.bed', required=True)
    group.add_argument('-p', '--pair', help='pair list from chromosomes: qry_id, ref_id', required=True)
    
    return group.parse_args()

# Longest crease sub sequence
def LIS(arr):
    min_num = [-1]
    for n in arr:
        k = bisect.bisect_left(min_num, n)
        if len(min_num) == k:
            min_num.append(n)
        else:
            min_num[k] = n
    return len(min_num)-1


def eval_synteny(ref_file, qry_file, pair_file):
    print("Loading bed files")
    ref_db = {}
    qry_db = {}
    with open(ref_file, 'r') as fin:
        for line in fin:
            data = line.strip().split()
            chrn = data[0]
            gid = data[3]
            if chrn not in ref_db:
                ref_db[chrn] = []
            ref_db[chrn].append(gid)
    
    with open(qry_file, 'r') as fin:
        for line in fin:
            data = line.strip().split()
            chrn = data[0]
            gid = data[3]
            if chrn not in qry_db:
                qry_db[chrn] = []
            qry_db[chrn].append(gid)
    
    print("Loading pairs")
    pair_list = []
    with open(pair_file, 'r') as fin:
        for line in fin:
            if line.strip() == '':
                continue
            pair_list.append(line.strip().split())
    
    print("Comparing")
    total_lis_values = 0
    total_gene_cnt = 0

    print("Result")
    for pair in pair_list:
        qchr, rchr = pair
        ref_idx_db = {ref_db[rchr][idx]: idx for idx in range(len(ref_db[rchr]))}
        tmp_order_list = [ref_idx_db[qry_db[qchr][idx]] 
                             if qry_db[qchr][idx] in ref_idx_db 
                             else -1 
                             for idx in range(len(qry_db[qchr]))]
        region_order_list = []
        for _ in tmp_order_list:
            if _ != -1:
                region_order_list.append(_)
        chr_gene_cnt = (len(ref_db[rchr])+len(qry_db[qchr]))/2.
        lis_value = max(LIS(region_order_list), LIS(region_order_list[::-1]))
        print("%s: %.4f%%"%(qchr, lis_value*100./chr_gene_cnt))
        total_lis_values += lis_value
        total_gene_cnt += chr_gene_cnt
    print("Total: %.4f%%"%(total_lis_values*100./total_gene_cnt))


if __name__ == "__main__":
    opts = get_opts()
    ref_file = opts.ref
    qry_file = opts.qry
    pair_file = opts.pair
    eval_synteny(ref_file, qry_file, pair_file)
