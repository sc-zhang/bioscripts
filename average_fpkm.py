#!/usr/bin/env python
import sys
import numpy as np


def average_fpkm(in_fpkm, out_avg):
	group_list = []
	last_smp = ''
	with open(in_fpkm, 'r') as fin:
		with open(out_avg, 'w') as fout:
			for line in fin:
				data = line.strip().split()
				if data[0] == 'gene_id':
					fout.write("gene_id")
					for i in range(1, len(data)):
						smp = data[i][:-1]
						if smp != last_smp:
							last_smp = smp
							group_list.append([])
							fout.write("\t%s"%smp)
						group_list[-1].append(i)
					fout.write("\n")
				else:
					fout.write("%s"%data[0])
					for idxs in group_list:
						vals = []
						for idx in idxs:
							vals.append(float(data[idx]))
						fout.write("\t%.2f"%np.average(vals))
					fout.write("\n")


if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Usage: python %s <in_fpkm> <out_avg>"%sys.argv[0])
	else:
		in_fpkm, out_avg = sys.argv[1:]
		average_fpkm(in_fpkm, out_avg)
