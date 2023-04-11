#!/usr/bin/env python
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import gzip
import argparse
import re
import time


def get_opts():
	group = argparse.ArgumentParser()
	group.add_argument("-b", "--base", help="Input vcf file as base, .gz supported", required=True)
	group.add_argument("-v", "--validation", help="Input vcf file as validation, .gz supported", required=True)
	group.add_argument("-r", "--repeat", help="Repeat regions file, gff format", default="")
	group.add_argument("-o", "--output", help="Output vcf file based on base vcf file, compressed with gzip", required=True)
	group.add_argument("-m", "--missing_rate", type=float, help="Missing rate threshold, percentage, default: 40", default=40)
	group.add_argument("-d", "--min_distance", type=int, help="Minimum distance between two snp sites, default: 0", default=0)
	return group.parse_args()


def read_vcf(in_vcf, method):
	if in_vcf[-3:].lower() == '.gz':
		f_in = gzip.open(in_vcf, 'rt')
	else:
		f_in = open(in_vcf, 'r')
	
	header = []
	vcf_infos = {}
	dp_db = {}
	for line in f_in:
		if line[0] == '#':
			header.append(line)
		else:
			data = line.strip().split()
			chrn = data[0]
			pos = int(data[1])
			ref = data[3]
			alt = data[4]
			filter = data[6]
			
			# filter indels
			if len(ref) > 1 or len(alt) > 1:
				continue
			if method.lower() == 'full':
				# filter LowQual marker
				if filter.lower() == 'lowqual':
					continue
			if chrn not in vcf_infos:
				vcf_infos[chrn] = {}
			if method.lower() == 'full':
				# calc missing rate
				m_c = 0
				for i in range(9, len(data)):
					if data[i].split(':')[0] == './.':
						m_c += 1
				m_r = m_c*1.0/(len(data)-8)
				vcf_infos[chrn][pos] = {'alt': alt, 'line': line, 'mr': m_r}
				dp = int(data[7].split('DP=')[1].split(';')[0])
				if dp not in dp_db:
					dp_db[dp] = 0
				dp_db[dp] += 1
			else:
				vcf_infos[chrn][pos] = {'alt': alt}
	f_in.close()
	if method.lower() == 'full':
		return header, vcf_infos, dp_db
	else:
		return vcf_infos


def read_gff(in_gff):
	regions_db = {}
	with open(in_gff, 'r') as fin:
		for line in fin:
			if line[0] == '#':
				continue
			data = line.strip().split()
			chrn = data[0]
			sr = int(data[3])
			er = int(data[4])
			if sr > er:
				tmp = er
				er = sr
				sr = tmp
			if chrn not in regions_db:
				regions_db[chrn] = []
			regions_db[chrn].append([sr, er])
	
	for chrn in regions_db:
		regions_db[chrn] = sorted(regions_db[chrn])

	return regions_db


def is_repeat(repeats, pos):
	s = 0
	e = len(repeats)-1
	while s<=e:
		mid = int((s+e)/2)
		if repeats[mid][0] > pos:
			e = mid-1
		elif repeats[mid][0] < pos:
			s = mid+1
		else:
			return True
	if repeats[e][1] >= pos:
		return True
	else:
		return False



def snp_filter(in_base, in_valid, in_rep, mr, md, out_file):
	print("\033[32m%s\033[0m Reading valid file"%(time.strftime('[%H:%M:%S]',time.localtime(time.time()))))
	valid_snps = read_vcf(in_valid, 'simple')

	print("\033[32m%s\033[0m Reading base file"%(time.strftime('[%H:%M:%S]',time.localtime(time.time()))))
	base_header, base_snps, base_dps = read_vcf(in_base, 'full')

	print("\033[32m%s\033[0m Reading repeat file"%(time.strftime('[%H:%M:%S]',time.localtime(time.time()))))
	if in_rep != '':
		repeat_db = read_gff(in_rep)
	else:
		repeat_db = {}
	'''
	dp_cnt = []
	for dp in sorted(base_dps):
		dp_cnt.append(base_dps[dp])
	dp_cnt_th = int(max(dp_cnt)*0.05)
	
	dp_min = -1
	for dp in sorted(base_dps):
		if base_dps[dp] > dp_cnt_th:
			if dp_min == -1:
				dp_min = dp
			dp_max = dp
	'''
	# Plot dist
	plt.figure(figsize=(10, 8), dpi=100)
	dp_x = []
	dp_y = []
	dp_total = []
	for dp in sorted(base_dps):
		plt.bar(x=dp, height=base_dps[dp], width=1, edgecolor='white', facecolor='blue', align='center', linewidth=0.01)
		dp_x.append(dp)
		dp_y.append(base_dps[dp])
		for i in range(0, base_dps[dp]):
			dp_total.append(dp_x)
	plt.plot(dp_x, dp_y, linewidth=0.05, linestyle='-', markersize=0, marker=',')

	dp_min = 5 #mean-std*1.96
	sum = np.sum(dp_y)
	top_sum = sum*0.95
	cnt = 0
	
	for dp in sorted(base_dps):
		if dp < dp_min:
			continue
		cnt += base_dps[dp]
		if cnt >= top_sum:
			dp_max = dp
			break
	plt.xlim(dp_x[0]-1, dp_max+1)
	plt.savefig('dist.pdf', filetype='pdf', bbox_inches='tight')	
	print("\033[32m%s\033[0m range=[%f, %f]"%(time.strftime('[%H:%M:%S]',time.localtime(time.time())), dp_min, dp_max))	

	print("\033[32m%s\033[0m Filtering and writing results"%(time.strftime('[%H:%M:%S]',time.localtime(time.time()))))
	if out_file[-3:].lower() != '.gz':
		out_file += '.gz'
	with gzip.open(out_file, 'wt') as fout:
		fout.write(''.join(base_header))
		for chrn in sorted(base_snps):
			if chrn not in valid_snps:
				continue
			last_pos = -1
			for pos in sorted(base_snps[chrn]):
				# filter repeat regions
				if chrn in repeat_db and is_repeat(repeat_db[chrn], pos):
					continue
				# filter base snps with valid snps
				if pos not in valid_snps[chrn]:
					continue
				if base_snps[chrn][pos]['alt'] != valid_snps[chrn][pos]['alt']:
					continue
				# filter missing rate
				if base_snps[chrn][pos]['mr'] > mr:
					continue
				
				# filter dp
				data = base_snps[chrn][pos]['line'].strip().split()
				dp = int(data[7].split('DP=')[1].split(';')[0])
				if dp < dp_min or dp > dp_max:
					continue
				if last_pos != -1 and pos - last_pos <= md:
					last_pos = pos
					continue
				else:
					last_pos = pos	
				fout.write(base_snps[chrn][pos]['line'])
	
	print("\033[32m%s\033[0m Finished"%(time.strftime('[%H:%M:%S]',time.localtime(time.time()))))
				

if __name__ == "__main__":
	opts = get_opts()
	in_base = opts.base
	in_valid = opts.validation
	in_rep = opts.repeat
	mr = opts.missing_rate/100.0
	md = opts.min_distance
	out_file = opts.output
	snp_filter(in_base, in_valid, in_rep, mr, md, out_file)
