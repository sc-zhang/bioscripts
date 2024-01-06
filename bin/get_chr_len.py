#!/usr/bin/env python
import sys
import gzip


def get_chr_len(in_fasta, out_file, is_chr_only):
	if is_chr_only != "T" and is_chr_only != "F":
		print("Error argument")
		exit(0)
	dict_len = {}
	if in_fasta[-3:].lower() == '.gz':
		f_in = gzip.open(in_fasta, 'rt')
	else:
		f_in = open(in_fasta, 'r')
	chrn = ''
	seq = ''
	for line in f_in:
		if line[0] == '>':
			if seq != '':
				dict_len[chrn] = len(seq)
			chrn = line.strip().split()[0][1:]
			seq = ''
		else:
			seq += line.strip()
	dict_len[chrn] = len(seq)
	f_in.close()
	
	with open(out_file, 'w') as f_out:
		chr_list = sorted(dict_len.keys())
		for chrn in chr_list:
			if is_chr_only == "T" and chrn[:3].lower() != 'chr':
				continue
			else:
				f_out.write(chrn+"\t"+str(dict_len[chrn])+"\n")


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Notice: script for calculating length of chromosomes in fasta file")
		print("Usage: python "+sys.argv[0]+" <fasta_file> <output_file> <T/F chr only>")
	else:
		f_fasta = sys.argv[1]
		f_out = sys.argv[2]
		is_chr_only = sys.argv[3]
		get_chr_len(f_fasta, f_out, is_chr_only)
