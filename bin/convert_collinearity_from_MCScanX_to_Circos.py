#!/usr/bin/env python
import sys


def get_col(in_col, in_gff, out_txt):
	id_db = {}
	with open(in_gff, 'r') as f_gff:
		for line in f_gff:
			data = line.strip().split()
			id_db[data[1]] = data[0][:2]+"Chr"+data[0][2:]+"\t"+data[2]+"\t"+data[3]+"\n"
	i = 0
	with open(in_col, 'r') as f_col:
		with open(out_txt, 'w') as f_out:
			for line in f_col:
				if line[0] == "#":
					continue
				data = line.strip().split('\t')
				id_1 = data[1]
				id_2 = data[2]
				if id_1 not in id_db or id_2 not in id_db:
					continue
				f_out.write("link"+str(i)+"\t"+id_db[id_1])
				f_out.write("link"+str(i)+"\t"+id_db[id_2])
				i += 1
			

if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Notice: script for converting collinearity file from MCScanX result to link file for Circos")
		print("Usage: python " + sys.argv[0] + " <collinearity_file> <gff_file> <out_file>")
	else:
		in_col = sys.argv[1]
		in_gff = sys.argv[2]
		out_txt = sys.argv[3]
		get_col(in_col, in_gff, out_txt)

