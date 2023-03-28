#!/usr/bin/env python
import sys


def get_tig_pos_of_chr(in_file):
	pos_db = {}
	with open(in_file, 'r') as f_in:
		for line in f_in:
			if line.strip() == '':
				continue
			data = line.strip().split()
			if data[4] == 'U':
				continue
			chrn = data[0]
			spos = int(data[1])
			epos = int(data[2])
			tig = data[5]
			direct = data[-1]
			pos_db[tig] = [chrn, spos, epos, direct]
	return pos_db


def convert_QTL_info(in_QTL, in_agp, out_QTL):
	tig_on_chr = get_tig_pos_of_chr(in_agp)
	with open(in_QTL, 'r') as fin:
		with open(out_QTL, 'w') as fout:
			for line in fin:
				data = line.strip().split()
				if data[0] == 'Pop':
					data.extend(['ActChr', 'Left_Pos', 'Right_Pos', 'Direct'])
				else:
					ltig, lpos = data[4].split('_')
					rtig, rpos = data[5].split('_')
					lpos = int(lpos)
					rpos = int(rpos)
					lchr, lsp, lep, ld = tig_on_chr[ltig]
					rchr, rsp, rep, rd = tig_on_chr[rtig]
					if lchr != rchr:
						print(data[1]+"\t"+lchr+"\t"+rchr)
						continue
					if ld == '-':
						lpos = lep-lpos+1
					else:
						lpos = lsp+lpos-1
					if rd == '-':
						rpos = rep-rpos+1
					else:
						rpos = rsp+rpos-1
					if lpos > rpos:
						tmp = lpos
						lpos = rpos
						rpos = tmp
						direct = "-"
					else:
						direct = "+"
					data.extend([lchr, str(lpos), str(rpos), direct])
				fout.write("%s\n"%'\t'.join(data))


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python "+sys.argv[0]+" <in_QTL> <in_agp> <out_QTL>")
	else:
		in_QTL, in_agp, out_QTL = sys.argv[1:]
		convert_QTL_info(in_QTL, in_agp, out_QTL)
