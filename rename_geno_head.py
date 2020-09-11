#!/usr/bin/env python
import sys


def rename_vcf_head(in_vcf, in_list, out_vcf):
	species_db = {}
	with open(in_list, 'r') as f_in:
		for line in f_in:
			data = line.strip().split()
			species_db[data[0]] = data[1]
	
	with open(in_vcf, 'r') as f_in:
		with open(out_vcf, 'w') as f_out:
			for line in f_in:
				if line[0] == "#":
					data = line.strip().split()
					for i in range(2, len(data)):
						if data[i] in species_db:
							data[i] = species_db[data[i]]
					f_out.write("\t".join(data)+"\n")
				else:
					f_out.write(line)


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Notice: modify head line in geno file with list")
		print("Usage: python "+sys.argv[0]+" <in_geno> <in_list> <out_geno>")
	else:
		in_vcf = sys.argv[1]
		in_list = sys.argv[2]
		out_vcf = sys.argv[3]
		rename_vcf_head(in_vcf, in_list, out_vcf)

