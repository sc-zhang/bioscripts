#!/usr/bin/env python
import sys, os


def bin_search(regions, v):
	s = 0
	e = len(regions) - 1
	while s <= e:
		m = (s+e)/2
		if regions[m][0] < v:
			s = m+1
		elif regions[m][0] > v:
			e = m-1
		else:
			return True
	if regions[e][0] > v:
		e -= 1
	if regions[e][0] <= v <= regions[e][1]:
		return True
	else:
		return False


def search(regions, v):
	for region in regions:
		if region[0] <= v <= region[1]:
			return True
	return False


def filter_snp(in_snp, result_dir, out_snp):
	if result_dir[0] != '/' or result_dir[0] != '.':
		result_dir = os.path.join(os.getcwd(), result_dir)
	files = os.listdir(result_dir)
	remove_db = []
	for fn in files:
		with open(result_dir+"/"+fn, 'r') as f_in:
			for line in f_in:
				data = line.strip().split(',')
				if len(data) < 8:
					continue
				if data[8] != 'nan' and data[8] != 'D':
					win = [int(data[1]),int(data[2])]
					if win not in remove_db:
						remove_db.append(win)
	remove_db = sorted(remove_db)
	
	with open(in_snp, 'r') as f_in:
		with open(out_snp, 'w') as f_out:
			for line in f_in:
				if line[0] == '#':
					f_out.write(line)
				else:
					pos = int(line.strip().split()[1])
					if search(remove_db, pos):
						continue
					else:
						f_out.write(line)


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Notice: script for removing snp position in regions of ABBA")
		print("Usage: python "+sys.argv[0]+" <in_snp> <result_dir> <out_snp>")
	else:
		in_snp = sys.argv[1]
		result_dir = sys.argv[2]
		out_snp = sys.argv[3]
		filter_snp(in_snp, result_dir, out_snp)

