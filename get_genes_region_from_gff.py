#!/usr/bin/env python
import sys


def get_genes_region_from_gff(gene_list, in_gff, out_bed):
	genes = []
	with open(gene_list, 'r') as f_in:
		for line in f_in:
			genes.append(line.strip())
	
	with open(in_gff, 'r') as f_in:
		with open(out_bed, 'w') as f_out:
			for line in f_in:
				if line[0] == '#' or line.strip() == '':
					continue
				data = line.strip().split()
				if data[2] != 'gene':
					continue
				id = data[8].split(';')[1].split('=')[1]
				if id in genes:
					f_out.write(data[0]+'\t'+data[3]+'\t'+data[4]+'\t'+id+'\n')


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python "+sys.argv[0]+" <gene_list> <in_gff> <out_bed>")
	else:
		proc, gene_list, in_gff, out_bed = sys.argv
		get_genes_region_from_gff(gene_list, in_gff, out_bed)

