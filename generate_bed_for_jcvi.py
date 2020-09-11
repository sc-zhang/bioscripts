#!/usr/bin/python
import sys
import re


def generate_bed_for_jcvi(in_gff, out_bed, chr_prefix):
	with open(in_gff, 'r') as f_in:
		with open(out_bed, 'w') as f_out:
			for line in f_in:
				if line[0] == "#":
					continue
				if line.startswith(chr_prefix) == False:
					continue
				data = line.strip().split()
				if data[2] != 'gene':
					continue
				chrn = data[0]
				#if chrn[-1].isdigit() == False:
				#	continue
				chr_pre = chrn[:len(chr_prefix)]
				chr_n = chrn[len(chr_prefix):]
				if len(chr_n) == 1:
					chrn = chr_pre+'0'+chr_n
				sp = int(data[3])
				ep = int(data[4])
				if sp == ep:
					continue
				if sp > ep:
					tmp = sp
					sp = ep
					ep = tmp
				region = "%d\t%d"%(sp, ep)
				direct = data[6]
				id = data[8].split(';')[1].split('=')[1]
				#f_out.write(chrn+'\t'+region+'\t'+id.replace('.', '_')+'\t0\t'+direct+'\n')
				f_out.write(chrn+'\t'+region+'\t'+id+'\t0\t'+direct+'\n')


if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Notice: this script is used to extract gene region from gff file for jvci module")
		print("Usage python "+sys.argv[0]+" <in_gff> <out_bed> [<chr_prefix>]")
	else:
		in_gff = sys.argv[1]
		out_bed = sys.argv[2]
		if len(sys.argv) < 4:
			chr_prefix = 'Chr'
		else:
			chr_prefix = sys.argv[3]
		generate_bed_for_jcvi(in_gff, out_bed, chr_prefix)

