#!/usr/bin/env python
import sys


def calc_gap_cnt(in_fa):
	fa_db = {}
	with open(in_fa, 'r') as fin:
		for line in fin:
			if line[0] == '>':
				id = line.strip()[1:]
				fa_db[id] = []
			else:
				fa_db[id].append(line.strip().lower())

	total_cnt = 0	
	for id in sorted(fa_db):
		last_base = ''
		cnt = 0
		seq = ''.join(fa_db[id])
		for i in range(0, len(seq)):
			if seq[i] == 'n' and last_base != 'n':
				cnt += 1
			last_base = seq[i]
		print("%s\t%d"%(id, cnt))
		total_cnt += cnt
	print("Total\t%d"%total_cnt)


if __name__ == "__main__":
	if len(sys.argv) < 2:
		print("Usage: python "+sys.argv[0]+" <in_fa>")
	else:
		in_fa = sys.argv[1]
		calc_gap_cnt(in_fa)
