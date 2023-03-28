#!/usr/bin/env python
import sys


def find_ovlp(in_gff3, out_bed):
	gff3_db = {}
	with open(in_gff3, 'r') as fin:
		for line in fin:
			if line[0] == '#' or line.strip() == '':
				continue
			data = line.strip().split()
			if data[2] != 'gene':
				continue
			chrn = data[0]
			sp = int(data[3])
			ep = int(data[4])
			dir = data[6]
			gn = data[8].split(";")[1].split("=")[1]
			if chrn not in gff3_db:
				gff3_db[chrn] = {}
			if dir not in gff3_db[chrn]:
				gff3_db[chrn][dir] = []
			gff3_db[chrn][dir].append([sp, ep, gn])
	
	with open(out_bed, 'w') as fout:
		for chrn in sorted(gff3_db):
			for dir in sorted(gff3_db[chrn]):
				tmp_list = []
				pos_list = sorted(gff3_db[chrn][dir])
				for i in range(0, len(pos_list)):
					if len(tmp_list) == 0:
						tmp_list.append(pos_list[i])
					else:
						s, e, gn = pos_list[i]
						if s <= last_e:
							tmp_list.append(pos_list[i])
						else:
							if len(tmp_list) > 1:
								for s, e, gn in tmp_list:
									fout.write("%s\t%d\t%d\t%s\t%s\n"%(chrn, s, e, dir, gn))
								fout.write("###\n")
							tmp_list = []
							tmp_list.append(pos_list[i])
					last_s, last_e, gn = pos_list[i]


if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Usage: python "+sys.argv[0]+" <in_gff3> <out_bed>")
	else:
		in_gff3, out_bed = sys.argv[1:]
		find_ovlp(in_gff3, out_bed)
