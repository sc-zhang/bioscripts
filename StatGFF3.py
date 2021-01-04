#!/usr/bin/env python3
import sys


def stat_gff3(in_gff3):
	anchored_db = {}
	unanchored_cnt = 0
	with open(in_gff3, 'r') as fin:
		for line in fin:
			if line.strip() == '' or line[0] == '#':
				continue
			data = line.strip().split()
			if data[2]!='gene':
				continue
			chrn = data[0]
			if chrn[:3]!='tig':
				allele = chrn[-1]
				chrn = chrn[:-1]
				if chrn not in anchored_db:
					anchored_db[chrn] = {}
				if allele not in anchored_db[chrn]:
					anchored_db[chrn][allele] = 0
				anchored_db[chrn][allele] += 1
			else:
				unanchored_cnt += 1
	
	for chrn in anchored_db:
		break

	print("Chromosome\t%s"%('\t'.join(sorted(anchored_db[chrn]))))
	anchored_cnt = 0
	
	for chrn in sorted(anchored_db):
		print("%s"%chrn, end='')
		for allele in sorted(anchored_db[chrn]):
			anchored_cnt += anchored_db[chrn][allele]
			print("\t%s"%("{:,}".format(anchored_db[chrn][allele])), end='')
		print("")
	print("Total no. of genes\t%s"%("{:,}".format(anchored_cnt+unanchored_cnt)))
	print("Unanchored genes\t%s"%("{:,}".format(unanchored_cnt)))
	print("Anchored genes\t%s"%("{:,}".format(anchored_cnt)))


if __name__ == "__main__":
	if len(sys.argv) < 2:
		print("Usage: python %s <in_gff3>"%sys.argv[0])
	else:
		in_gff3 = sys.argv[1]
		stat_gff3(in_gff3)
