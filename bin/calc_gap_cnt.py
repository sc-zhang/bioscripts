#!/usr/bin/env python
import sys


def calc_gap_cnt(in_fa):
	total_cnt = 0
	with open(in_fa, 'r') as fin:
		id = ''
		gap_cnt = 0
		for line in fin:
			if line[0] == '>':
				if id != '':
					print("%s\t%d"%(id, gap_cnt))
					total_cnt += gap_cnt
				id = line.strip()[1:]
				gap_cnt = 0
				last_base =''
			else:
				for i in range(len(line.strip())):
					if line[i].lower() == 'n' and last_base != 'n':
						gap_cnt += 1
					last_base = line[i].lower()
		print("%s\t%d"%(id, gap_cnt))
		total_cnt += gap_cnt
	print("Total\t%d"%total_cnt)


if __name__ == "__main__":
	if len(sys.argv) < 2:
		print("Usage: python "+sys.argv[0]+" <in_fa>")
	else:
		in_fa = sys.argv[1]
		calc_gap_cnt(in_fa)
