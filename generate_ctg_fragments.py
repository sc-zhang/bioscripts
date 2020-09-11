#!/usr/bin/python
import sys


def substr(source, n):
	l_s = len(source) / n
	strs = []
	for i in range(0, n):
		if i == n-1:
			strs.append(source[l_s*i:])
		else:
			strs.append(source[l_s*i: l_s*(i+1)])
	return strs


def split_ctg(in_ctg, out_ctg, frg_count):
	with open(in_ctg, 'r') as f_in:
		with open(out_ctg, 'w') as f_out:
			seq = ''
			id = ''
			for line in f_in:
				if line[0] == '>':
					if seq != '':
						seqs = substr(seq, frg_count)
						for i in range(0, len(seqs)):
							f_out.write(id+"_fragment"+str(i)+"\n"+seqs[i]+"\n")
					seq = ''
					id = line.strip()
				else:
					seq += line.strip()
			seqs = substr(seq, frg_count)
			for i in range(0, len(seqs)):
				f_out.write(id+"_fragment"+str(i)+"\n"+seqs[i]+"\n")


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Notice: split sequence into serval fragments")
		print("Usage: python "+sys.argv[0]+" <in_ctg_fasta> <out_ctg_fasta> <fragment_counts>")
	else:
		in_ctg = sys.argv[1]
		out_ctg = sys.argv[2]
		frg_count = int(sys.argv[3])
		split_ctg(in_ctg, out_ctg, frg_count)

