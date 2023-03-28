#!/usr/bin/python
import sys, os


def split_fasta_by_chr(fasta_file, out_folder):
	if os.path.exists(out_folder) == False:
		os.mkdir(out_folder)
	seq_db = {}
	with open(fasta_file, 'r') as f_fasta:
		seq = ''
		seq_id = ''
		for line in f_fasta:
			if line[0] == ">":
				if seq != '':
					seq_db[seq_id] = seq
				seq_id = line.strip()
				seq = ''
			else:
				seq += line.strip()
		seq_db[seq_id] = seq
	
	for seq_id in seq_db:
		if seq_id[:4].lower() != '>chr':
			continue
		f_out = open(out_folder+"/"+seq_id[1:]+".fasta", 'w')
		f_out.write(seq_id+"\n"+seq_db[seq_id])
		f_out.close()				


if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Notice: script for spliting fasta into serval files contain single chromosome")
		print("Usage: python " + sys.argv[0] + " <in_fasta> <out_dir>")
	else:
		in_fasta = sys.argv[1]
		out_dir = sys.argv[2]
		split_fasta_by_chr(in_fasta, out_dir)
