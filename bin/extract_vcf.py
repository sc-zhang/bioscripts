#!/usr/bin/env python
import sys


def bin_search(bed_list, pos):
	s = 0
	e = len(bed_list)-1
	while s <= e:
		mid = (s+e)/2
		if bed_list[mid][0] < pos:
			s = mid+1
		elif bed_list[mid][0] > pos:
			e = mid-1
		else:
			return True
	if bed_list[e][0] <= pos and bed_list[e][1] >= pos:
		return True
	else:
		return False


def extract_vcf(in_vcf, in_bed, out_vcf):
	bed_db = {}
	with open(in_bed, 'r') as f_in:
		for line in f_in:
			if line.strip() == "":
				continue
			data = line.strip().split()
			chrn = data[0]
			sr = int(data[1])
			er = int(data[2])
			if chrn not in bed_db:
				bed_db[chrn] = []
			bed_db[chrn].append([sr, er])
	with open(in_vcf, 'r') as f_in:
		with open(out_vcf, 'w') as f_out:
			for line in f_in:
				if line.strip() == "":
					continue
				if line[0] == '#':
					f_out.write(line)
				else:
					data = line.strip().split()
					chrn = data[0]
					pos = int(data[1])
					if bin_search(bed_db[chrn], pos):
						f_out.write(line)


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python "+sys.argv[0]+" <in_vcf> <in_bed> <out_vcf>")
	else:
		proc, in_vcf, in_bed, out_vcf = sys.argv
		extract_vcf(in_vcf, in_bed, out_vcf)

