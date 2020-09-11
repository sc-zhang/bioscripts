#!/usr/bin/env python
import sys


def reverse_comp(seq):
	rev_seq = ''
	for i in range(len(seq)-1, -1, -1):
		if seq[i].lower() == 'a':
			rev_seq += 'T' 
		elif seq[i].lower() == 't':
			rev_seq += 'A' 
		elif seq[i].lower() == 'g':
			rev_seq += 'C' 
		elif seq[i].lower() == 'c':
			rev_seq += 'G' 
		else:
			rev_seq += seq[i]
	return rev_seq


def extract_gene_seq(in_fa, in_list, out_fa):
	with open(in_fa, 'r') as f_in:
		seq = ''
		for line in f_in:
			if line[0] == '>':
				continue
			else:
				seq += line.strip()
	
	ext_db = {}
	with open(in_list) as f_in:
		for line in f_in:
			data = line.strip().split()
			gs = int(data[0])
			ge = int(data[1])
			gn = data[2]
			direct = data[3]
			if gn not in ext_db:
				ext_db[gn] = []
			ext_db[gn].append([gs, ge, direct])

	with open(out_fa, 'w') as f_out:
		for gn in ext_db:
			for i in range(0, len(ext_db[gn])):
				gs, ge, direct = ext_db[gn][i]
				if direct == '+':
					gene_seq = seq[gs: ge]
				else:
					gene_seq = reverse_comp(seq[gs: ge])
				if i == 0:
					gene = gn
				else:
					gene = gn+"_"+str(i)
				f_out.write(">%s\n%s\n"%(gene, gene_seq))


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python "+sys.argv[0]+" <in_fasta> <in_list> <out_fasta>")
	else:
		in_fa, in_list, out_fa = sys.argv[1:]
		extract_gene_seq(in_fa, in_list, out_fa)
