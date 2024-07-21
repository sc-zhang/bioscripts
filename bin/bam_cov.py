#!/usr/bin/env python3
import pysam
import multiprocessing
import argparse


def get_opts():
    group = argparse.ArgumentParser()
    group.add_argument("-b", "--bam", help="Input bam file, must be indexed", required=True)
    group.add_argument("-o", "--output", help="Output statistic", required=True)
    group.add_argument("-t", "--threads", help="Threads, default=10", type=int, default=10)
    return group.parse_args()


def sub_cov(in_bam, chrn, chrl):
    bins = [0 for _ in range(chrl)]
    with pysam.AlignmentFile(in_bam, 'rb') as fin:
        for read in fin.fetch(chrn):
            for pos in read.get_reference_positions():
                bins[pos] = 1
    return chrn, chrl, sum(bins)


def main():
    opts = get_opts()
    in_bam = opts.bam
    out_stat = opts.output
    threads = opts.threads
    chr_name = []
    chr_len = []
    with pysam.AlignmentFile(in_bam, 'rb') as fin:
        chr_name = fin.references
        chr_len = fin.lengths
    
    res = []
    pool = multiprocessing.Pool(processes=threads)
    for _ in range(len(chr_name)):
        r = pool.apply_async(sub_cov, (in_bam, chr_name[_], chr_len[_], ))
        res.append(r)
    pool.close()
    pool.join()

    with open(out_stat, 'w') as fout:
        total_covl = 0
        total_chrl = 0
        for r in res:
            chrn, chrl, covl = r.get()
            fout.write("%s\t%d\t%d\t%f\n"%(chrn, chrl, covl, covl*1./chrl))
            total_covl += covl
            total_chrl += chrl
        fout.write("Total\t%s\t%d\t%f\n"%(total_chrl, total_covl, total_covl*1./total_chrl))


if __name__ == "__main__":
    main()
