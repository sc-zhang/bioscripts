#!/usr/bin/env python
import sys, os


def split_fasta_by_count(in_fa, is_seq, cnt, out_dir):
	if not os.path.isdir(out_dir):
		os.mkdir(out_dir)
	
	seq_db = {}
	id = ''
	seq = ''
	with open(in_fa, 'r') as f_in:
		for line in f_in:
			if line[0] == '>':
				if seq != '':
					seq_db[id] = seq
				id = line.strip()[1:]
				seq = ''
			else:
				seq += line
		seq_db[id] = seq
	
	total_seq_cnt = len(seq_db)
	tmp_cnt = int(round(total_seq_cnt*1.0/cnt+0.5))
	if is_seq:
		file_cnt = tmp_cnt
		seq_cnt = cnt
	else:
		file_cnt = cnt
		seq_cnt = tmp_cnt
	fn = in_fa.replace('.fasta', '').replace('.fa', '')
	id_list = seq_db.keys()
	for i in range(0, file_cnt):
		with open(os.path.join(out_dir, fn+"_"+str(i)+".fa"), 'w') as f_out:
			for j in range(0, seq_cnt):
				index = i*seq_cnt+j
				if index < len(id_list):
					f_out.write(">%s\n%s"%(id_list[index], seq_db[id_list[index]]))


if __name__ == "__main__":
	if len(sys.argv) < 5:
		print("Usage: python "+sys.argv[0]+" <in_fasta> <S/F> <count> <out_dir>")
	else:
		in_fa = sys.argv[1]
		if sys.argv[2].lower() == 's':
			is_seq = True
		else:
			is_seq = False
		cnt = int(sys.argv[3])
		out_dir = sys.argv[4]
		split_fasta_by_count(in_fa, is_seq, cnt, out_dir)

