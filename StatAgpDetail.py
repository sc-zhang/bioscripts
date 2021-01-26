#!/usr/bin/env python3
import sys


def stat_agp(in_agp, out_csv):
	asm_db = {}
	unanc_cnt = 0
	unanc_len = 0
	anc_cnt = 0
	anc_len = 0
	gap_cnt = 0
	gap_len = 0
	with open(in_agp, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			if data[4] == 'U':
				gap_cnt += 1
				gap_len += 100
			else:
				chrn = data[0]
				if chrn[:3] != 'tig':
					allele = chrn[-1]
					chrn = chrn[:-1]
					if chrn not in asm_db:
						asm_db[chrn] = {}
					if allele not in asm_db[chrn]:
						asm_db[chrn][allele] = {'cnt': 0, 'len': 0}
					asm_db[chrn][allele]['cnt'] += 1
					asm_db[chrn][allele]['len'] = int(data[2])
					anc_cnt += 1
					anc_len += int(data[7])
				else:
					unanc_cnt += 1
					unanc_len += int(data[2])
	
	for chrn in asm_db:
		break
	with open(out_csv, 'w') as fout:
		fout.write(",%s\n"%(',,'.join(sorted(asm_db[chrn]))))
		for chrn in sorted(asm_db):
			info = [chrn]
			for allele in sorted(asm_db[chrn]):
				info.append("\"%s\""%("{:,}".format(asm_db[chrn][allele]['cnt'])))
				info.append("\"%s\""%("{:,}".format(asm_db[chrn][allele]['len'])))
			fout.write("%s\n"%(','.join(info)))
		fout.write("Anchored contigs,\"%s\",\"%s\"\n"%("{:,}".format(anc_cnt), "{:,}".format(anc_len/1e6)))
		fout.write("Unanchored contigs,\"%s\",\"%s\"\n"%("{:,}".format(unanc_cnt), "{:,}".format(unanc_len/1e6)))
		fout.write("Gaps,\"%s\",\"%s\"\n"%("{:,}".format(gap_cnt), "{:,}".format(gap_len/1e6)))
		fout.write("Total,\"%s\",\"%s\"\n"%("{:,}".format(anc_cnt+unanc_cnt), "{:,}".format((anc_len+unanc_len)/1e6)))
	

if __name__ == "__main__":
	if len(sys.argv) < 2:
		print("Usage: python %s <in_agp> <out_csv>"%sys.argv[0])
	else:
		in_agp, out_csv = sys.argv[1:]
		stat_agp(in_agp, out_csv)
