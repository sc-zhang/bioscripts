#!/usr/bin/env python
import sys


def extract_gene_from_gff(in_list, in_gff, out_bed):
	id_list = []
	with open(in_list, 'r') as f_in:
		for line in f_in:
			id_list.append(line.strip())
	
	with open(in_gff, 'r') as f_in:
		with open(out_bed, 'w') as f_out:
			for line in f_in:
				if line[0] == '#' or line.strip() == '':
					continue
				data = line.strip().split()
				if data[2] != 'gene':
					continue
				id = data[8].split(';')[0].split("=")[1]
				if id in id_list:
					f_out.write(data[0]+'\t'+data[3]+'\t'+data[4]+'\t'+id+'\n')


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Notice: this script is used to extract region of gene in list file from gff file")
		print("Usage: python "+sys.argv[0]+" <in_list> <in_gff> <out_bed>")
	else:
		in_list = sys.argv[1]
		in_gff = sys.argv[2]
		out_bed = sys.argv[3]
		extract_gene_from_gff(in_list, in_gff, out_bed)

