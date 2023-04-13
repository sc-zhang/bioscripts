#!/usr/bin/env python
import sys


def filter_cds(in_cds, out_cds):
	print("Loading cds")
	cds_db = {}
	with open(in_cds, 'r') as fin:
		for line in fin:
			if line[0] == '>':
				id = line.strip().split()[0][1:]
				cds_db[id] = []
			else:
				cds_db[id].append(line.strip().upper())
	
	for id in cds_db:
		cds_db[id] = ''.join(cds_db[id])
	print("Filtering cds")
	start_codon = set(["ATG"])
	stop_codon = set(["TAG", "TAA", "TGA"])

	with open(out_cds, 'w') as fout:
		for id in sorted(cds_db):
			cds_len = len(cds_db[id])
			cds_start = cds_db[id][:3]
			cds_stop = cds_db[id][-3:]
			if (cds_len%3 != 0) or (cds_start not in start_codon) or (cds_stop not in stop_codon):
				fout.write(">%s\n%s\n"%(id, cds_db[id]))
	
	print("Finished")
		

if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Usage: python %s <in_cds> <out_cds>"%sys.argv[0])
	else:
		in_cds, out_cds = sys.argv[1:]
		filter_cds(in_cds, out_cds)

