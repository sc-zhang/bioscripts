#!/usr/bin/env python
import sys


def merge_regions(dict_regions):
	coverage_db = {}
	for chrn in dict_regions:
		if len(dict_regions[chrn]) == 0:
			continue
		coverage_db[chrn] = []
		dict_regions[chrn] = sorted(dict_regions[chrn])
		last_s = 0
		last_e = 0
		temp_regions = []
		for i in range(0, len(dict_regions[chrn])):
			s = dict_regions[chrn][i][0]
			e = dict_regions[chrn][i][1]
			if i == 0:
				temp_regions.append(s)
				last_s = s
				last_e = e
				continue
			if s <= last_e:
				if e > last_e:
					last_e = e
			else:
				temp_regions.append(last_e)
				temp_regions.append(s)
				last_s = s
				last_e = e
		temp_regions.append(last_e)
		for i in range(0, len(temp_regions), 2):
			coverage_db[chrn].append([temp_regions[i],temp_regions[i+1]])
	return coverage_db


def filter_coverage(in_file, blast_result, outfile, threshold):
	chr_db = {}
	with open(in_file, 'r') as f_in:
		id = ''
		seq = ''
		for line in f_in:
			if line[0] == ">":
				if seq != '':
					chr_db[id] = seq
				id = line.strip()[1:]
				seq = ''
			else:
				seq += line.strip()
			chr_db[id] = seq
	

	regions = {}
	for chrn in chr_db:
		regions[chrn] = []
	
	with open(blast_result, 'r') as f_in:
		for line in f_in:
			if line[0] == '#':
				continue
			data = line.strip().split()
			chrn = data[0]
			if chrn not in regions:
				continue
			value = float(data[2])
			if value < threshold:
				continue
			sr = int(data[6])
			er = int(data[7])
			if sr > er:
				temp = sr
				sr = er
				er = temp
			regions[chrn].append([sr, er])
	
	coverages = merge_regions(regions)
	
	chr_list = sorted(coverages.keys())
	with open(out_file, 'w') as f_out:
		for chrn in chr_list:
			ss = 1#coverages[chrn][0][0]
			ee = len(chr_db[chrn])#coverages[chrn][len(coverages[chrn])-1][1]
			for i in range(0, len(coverages[chrn])):
				s = coverages[chrn][i][0]
				e = coverages[chrn][i][1]
				if ss != s:
					f_out.write('>'+chrn+'['+str(ss)+':'+str(s-1)+']\n')
					f_out.write(chr_db[chrn][ss-1:s-1]+'\n')
				ss = e+1
			if ss <= ee:
				f_out.write('>'+chrn+'['+str(ss)+':'+str(ee-1)+']\n')
				f_out.write(chr_db[chrn][ss-1: ee]+'\n')


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python "+sys.argv[0]+" <in_fasta> <blast_result> <out_file> [<threshold>]")
	else:
		in_file = sys.argv[1]
		blast_result = sys.argv[2]
		out_file = sys.argv[3]
		if len(sys.argv) == 5:
			threshold = float(sys.argv[4])
		else:
			threshold = 0.0
		filter_coverage(in_file, blast_result, out_file, threshold)

