#!/usr/bin/python
import sys, os


def split_fasta_by_chr(tbl_file, out_folder):
	if os.path.exists(out_folder) == False:
		os.mkdir(out_folder)
	seq_db = {}
	with open(tbl_file, 'r') as f_tbl:
		seq = ''
		seq_id = ''
		for line in f_tbl:
			if line[0] == ">":
				if seq != '':
					seq_db[seq_id] = seq
				seq_id = line.strip()
				seq = ''
			else:
				seq += line
		seq_db[seq_id] = seq
	
	for seq_id in seq_db:
		f_out = open(out_folder+"/"+seq_id.split()[-1]+".tbl", 'w')
		f_out.write(seq_id+"\n"+seq_db[seq_id])
		f_out.close()				


if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Notice: script for spliting tbl into serval files contain single chromosome")
		print("Usage: python " + sys.argv[0] + " <in_tbl> <out_dir>")
	else:
		in_tbl = sys.argv[1]
		out_dir = sys.argv[2]
		split_fasta_by_chr(in_tbl, out_dir)
