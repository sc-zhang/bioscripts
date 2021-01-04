#!/usr/bin/env python3
import sys


def stat_agp(in_agp):
	asm_db = {}
	total_tig = 0
	unchor_tig = 0
	unchor_tig_size = 0
	asm_size = 0
	with open(in_agp, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			if data[4] == 'U':
				continue
			chrn = data[0]
			total_tig += 1
			ep = int(data[2])
			if chrn[:3] != 'tig':
				allele = chrn[-1]
				chrn = chrn[:-1]
				if chrn not in asm_db:
					asm_db[chrn] = {}
				if allele not in asm_db[chrn]:
					asm_db[chrn][allele] = ep
				asm_db[chrn][allele] = ep
			else:
				unchor_tig += 1
				unchor_tig_size += ep
		for chrn in asm_db:
			for allele in asm_db[chrn]:
				asm_size += asm_db[chrn][allele]
	print("\t%s"%('\t'.join(sorted(asm_db[chrn]))))

	for chrn in sorted(asm_db):
		print("%s"%chrn, end='')
		for allele in sorted(asm_db[chrn]):
			print("\t%s"%("{:,}".format(asm_db[chrn][allele])), end='')
		print("")
	print("No. of unanchored contigs\t%s"%("{:,}".format(unchor_tig)))
	print("Unanchored sequences (Mb)\t%s"%("{:,}".format(unchor_tig_size*1.0/1e6)))
	print("Total no. of contigs\t%s"%("{:,}".format(total_tig)))
	print("Total assembled size (Mb)\t%s"%("{:,}".format(asm_size*1.0/1e6)))


if __name__ == "__main__":
	if len(sys.argv) < 2:
		print("Usage: python %s <in_agp>"%sys.argv[0])
	else:
		in_agp = sys.argv[1]
		stat_agp(in_agp)
