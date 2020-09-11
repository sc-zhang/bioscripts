#!/usr/bin/python
import sys


def get_seq(in_fasta, in_bed, out_fasta):
	seq_db = {}
	print("Loading fasta")
	with open(in_fasta, 'r') as f_in:
		id = ''
		seq = ''
		for line in f_in:
			if line[0] == '>':
				if seq != '':
					seq_db[id] = seq
				id = line.strip().split()[0][1:]
				seq = ''
			else:
				seq += line.strip()
		seq_db[id] = seq

	print("%d sequences loaded"%len(seq_db))
	
	print("Writing fasta")
	cnt = 0
	with open(in_bed, 'r') as f_in:
		with open(out_fasta, 'w') as f_out:
			for line in f_in:
				data = line.strip().split()
				chrn = data[0]
				s = int(data[1])
				e = int(data[2])
				if s>e:
					s, e = e, s
				gene = data[3]
				f_out.write(">"+gene+"\n"+seq_db[chrn][s:e+1]+"\n")
				cnt += 1
	print("%d sequences wrote"%cnt)
	print("Success")


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Notice: extract sequences with bed file")
		print("Usage: python "+sys.argv[0]+" <in_fasta> <in_bed> <out_fasta>")
	else:
		in_fasta = sys.argv[1]
		in_bed = sys.argv[2]
		out_fasta = sys.argv[3]
		get_seq(in_fasta, in_bed, out_fasta)

