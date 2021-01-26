#!/usr/bin/env python
import sys
import re


def convert_simple_for_circos(in_simple, in_gff3_files, out_link):
	print("Loading gff3")
	id_db = {}
	for in_gff3 in in_gff3_files.split(','):
		with open(in_gff3, 'r') as fin:
			for line in fin:
				if line.strip() == '' or line[0] == '#':
					continue
				data = line.strip().split()
				if data[2]!='gene':
					continue
				chrn = data[0]
				sp = int(data[3])
				ep = int(data[4])
				if 'Name' in data[8]:
					id = re.findall(r'Name=(.*)', data[8])[0].split(';')[0]
				else:
					id = re.findall(r'ID=(.*)', data[8])[0].split(';')[0]
				id_db[id] = [chrn, sp, ep]
	
	print("Loading and writing link")
	with open(in_simple, 'r') as fin:
		with open(out_link, 'w') as fout:
			for line in fin:
				data = line.strip().split()
				achrn = id_db[data[0]][0]
				asp = min(id_db[data[0]][1], id_db[data[1]][1])
				aep = max(id_db[data[0]][2], id_db[data[1]][2])
				bchrn = id_db[data[2]][0]
				bsp = min(id_db[data[2]][1], id_db[data[3]][1])
				bep = max(id_db[data[2]][2], id_db[data[3]][2])
				fout.write("%s\t%d\t%d\t%s\t%d\t%d\n"%(achrn, asp, aep, bchrn, bsp, bep))
	
	print("Finished")				


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python %s <in_simple> <in_gff3_files> <out_link>"%sys.argv[0])
	else:
		in_simple, in_gff3_files, out_link = sys.argv[1:]
		convert_simple_for_circos(in_simple, in_gff3_files, out_link)
