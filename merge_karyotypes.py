#!/usr/bin/python
import sys


def merge_karyotype(ARGVS):
	f_cl = ARGVS[0]
	f_out = ARGVS[1]

	f_L = []
	with open(f_cl, 'r') as f_in:
		for line in f_in:
			data = line.strip().split()
			for fn in data:
				if fn != '':
					f_L.append(fn)
	
	dict_chr_len = {}
	color_map = {}
	for fn in f_L:
		with open(fn, 'r') as f_in:
			fn_pre = fn.split('.')[0]
			for line in f_in:
				data = line.strip().split()
				chrn = "Chr"+data[0][3:].zfill(2)
				ID = fn_pre+"."+chrn
				info = fn_pre + line.strip()
				dict_chr_len[ID] = info
				color_map[ID] = data[0].lower()
	with open(f_out, 'w') as f_o:
		chr_key = sorted(dict_chr_len.keys())
		for chrn in chr_key:
			data = dict_chr_len[chrn].split()
			f_o.write("chr\t-\t"+data[0]+"\t"+chrn+"\t0\t"+data[1]+"\t"+color_map[chrn]+"\n")


if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Notice: script for merging karyotypes of different species working for Circos")
		print("Usage: python "+sys.argv[0]+" <chr_len_files> <outputfile>")
		print("    chr_len_files : a file contain filenames of chrlenfiles of species,")
		print("                    filenames should be start with short name of species")
		print("                    and need a dot after short name")
	else:
		ARGVS = sys.argv[1:]
		merge_karyotype(ARGVS)

