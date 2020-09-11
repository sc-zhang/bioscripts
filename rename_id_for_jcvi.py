#!/usr/bin/python
import sys


def rename_id(in_fasta, out_fasta):
	with open(in_fasta, 'r') as f_in:
		with open(out_fasta, 'w') as f_out:
			for line in f_in:
				if line[0] == '>':
					data = line.strip().split()[0]
					id = data.split('.')
					if len(id) > 1:
						id = id[0]+'_'+id[1]
					else:
						id = id[0]
					f_out.write(id+'\n')
				else:
					f_out.write(line)


if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Notice: script for rename id in fasta file for jvci module")
		print("Usage: python "+sys.argv[0]+" <in_fasta> <out_fasta>")
	else:
		in_fasta = sys.argv[1]
		out_fasta = sys.argv[2]
		rename_id(in_fasta, out_fasta)

