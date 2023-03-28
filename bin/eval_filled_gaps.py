#!/usr/bin/env python
import sys


def search_gaps(seq):
	gaps_db = []
	cnt_n = 0
	for i in range(0, len(seq)):
		if seq[i].lower() == 'n':
			if cnt_n == 0:
				s = i
			cnt_n += 1
		else:
			if cnt_n != 0:
				e = i
				gaps_db.append([s, e-1])
				cnt_n = 0
	cnt_n = 0
	for region in gaps_db:
		cnt_n += region[1]-region[0]+1
	return gaps_db, cnt_n


def calc_gaps(seq):
	cnt_n = 0
	for i in range(0, len(seq)):
		if seq[i].lower() == 'n':
			cnt_n += 1
	return cnt_n


def make_seq_db(in_fasta):
	seq_db = {}
	with open(in_fasta, 'r') as f_in:
		id = ''
		seq = ''
		for line in f_in:
			if line[0] == ">":
				if seq != '':
					seq_db[id] = seq
				id = line.strip()
				seq = ''
			else:
				seq += line.strip()
		seq_db[id] = seq
	return seq_db


def eval_filled_gaps(ref_fasta, query_fasta, result_file):
	print("Reading reference fasta")
	ref_seq_db = make_seq_db(ref_fasta)

	print("Reading query fasta")
	query_seq_db = make_seq_db(query_fasta)

	print("Evaluating")
	with open(result_file, 'w') as f_out:
		for id in ref_seq_db:
			ref_gaps_db, ref_gaps_cnt = search_gaps(ref_seq_db[id])
			query_gaps_cnt = calc_gaps(query_seq_db[id])
			f_out.write(id[1:]+"\n")
			if ref_gaps_cnt != 0:
				f_out.write("Filled %0.2f%%\n"%((ref_gaps_cnt-query_gaps_cnt)*1.0/ref_gaps_cnt*100.0))
			else:
				f_out.write("No gaps\n")
			if len(ref_gaps_db) != 0:
				for region in ref_gaps_db:
					s = region[0]
					e = region[1]
					f_out.write("Region %d-%d:\n%s\n"%(s, e, query_seq_db[id][s:e+1]))
			f_out.write("\n")
	print("Success")


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Notice: this script is used to evaulate status that gaps been filled")
		print("Usage: python "+sys.argv[0]+" <ref_fasta> <query_fasta> <result_file>")
	else:
		ref_fasta = sys.argv[1]
		query_fasta = sys.argv[2]
		result_file = sys.argv[3]
		eval_filled_gaps(ref_fasta, query_fasta, result_file)

