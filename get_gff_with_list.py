#!/usr/bin/env python
import sys
import re


def get_gff_with_list(in_gff, in_list, out_gff):
	print("Loading list")
	use_id = {}
	with open(in_list, 'r') as fin:
		for line in fin:
			use_id[line.strip()] = 1
	
	print("Filter gff3")
	with open(in_gff, 'r') as fin:
		with open(out_gff, 'w') as fout:
			fout.write("#gff-version 3\n")
			is_write = False
			for line in fin:
				if line.strip() == '' or line[0] == '#':
					continue
				data = line.strip().split()
				if data[2] == 'gene':
					if "Name" in data[8]:
						regexp = r'Name=(.*)'
					else:
						regexp = r'ID=(.*)'
					id = re.findall(regexp, data[8])[0].split(';')[0]
					if id in use_id:
						is_write = True
					else:
						is_write = False
				if is_write:
					fout.write(line)
	
	print("Finished")


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python %s <in_gff> <in_list> <out_gff>"%sys.argv[0])
	else:
		in_gff, in_list, out_gff = sys.argv[1:]
		get_gff_with_list(in_gff, in_list, out_gff)
