#!/usr/bin/env python
import sys, operator


def reverse_region(ori_region, chr_len_db):
	new_region = {}
	for chrn in ori_region:
		if chrn not in new_region:
			new_region[chrn] = [] 
		temp_region = []
		temp_region.append(0)
		for region in ori_region[chrn]:
			temp_region.append(region[0]-1)
			temp_region.append(region[1]+1)
		temp_region.append(chr_len_db[chrn]-1)
		i = 0
		while i < len(temp_region):
			new_region[chrn].append([temp_region[i], temp_region[i+1]])
			i += 2
	return new_region


def update_region(ori_region, region):
	temp_region = []
	for chr_region in ori_region:
		if chr_region[0] <= region[0] <= chr_region[1]:
			if chr_region not in temp_region:
				temp_region.append(chr_region)
		if chr_region[0] <= region[1] <= chr_region[1]:
			if chr_region not in temp_region:
				temp_region.append(chr_region)
	for chr_region in temp_region:
		if chr_region in ori_region:
			ori_region.remove(chr_region)
		if chr_region[0] < region[0] <= chr_region[1]:
			ori_region.append([chr_region[0], region[0]-1])
		if chr_region[0] <= region[1] < chr_region[1]:
			ori_region.append([region[1]+1, chr_region[1]])
	return ori_region


def remove_region(blast_results, chr_len_file, out_bed):
	chr_len_db = {}
	with open(chr_len_file, 'r') as f_in:
		for line in f_in:
			if line.strip() != '':
				data = line.strip().split()
				chr_len_db[data[0]] = int(data[1])
	out_region = {}
	is_first = True
	for blast_file in blast_results:
		if is_first:
			with open(blast_file, 'r') as f_in:
				blast_region = {}
				for line in f_in:
					data = line.strip().split()
					chrn = data[1]
					s = int(data[8])
					e = int(data[9])
					if s > e:
						temp = s
						s = e
						e = temp
					if chrn not in out_region:
						blast_region[chrn] = []
						blast_region[chrn].append([s, e])
					temp_list = []
					for i in range(0, len(blast_region[chrn])):
						if s <= blast_region[chrn][i][0] and e >= blast_region[chrn][i][0]:
							blast_region[chrn][i][0] = s
						if e >= blast_region[chrn][i][1] and s <= blast_region[chrn][i][1]:
							blast_region[chrn][i][1] = e
						if s > blast_region[chrn][i][1] or e < blast_region[chrn][i][0]:
							temp_list.append([s, e])
					for region in temp_list:
						blast_region[chrn].append(region)
					blast_region[chrn] = sorted(blast_region[chrn])
				for chrn in blast_region:
					temp_region = []
					last_e = 0
					for i in range(0, len(blast_region[chrn])):
						s = blast_region[chrn][i][0]
						e = blast_region[chrn][i][1]
						if i == 0:
							temp_region.append(s)
							last_e = e
						else:
							if last_e < s:
								temp_region.append(e)
								temp_region.append(s)
								last_e = e
							else:
								if last_e < e:
									last_e = e
					temp_region.append(e)
					i = 0
					blast_region[chrn] = []
					while i < len(temp_region):
						blast_region[chrn].append([temp_region[i], temp_region[i+1]])
						i += 2
					out_region = reverse_region(blast_region, chr_len_db)
			is_first = False
		else:
			with open(blast_file, 'r') as f_in:
				for line in f_in:
					data = line.strip().split()
					chrn = data[1]
					s = int(data[8])
					e = int(data[9])
					if s > e:
						temp = s
						s = e
						e = temp
					out_region[chrn] = update_region(out_region[chrn], [s, e])
	
	with open(out_bed, 'w') as f_out:
		for chrn in out_region:
			out_region[chrn] = sorted(out_region[chrn])
			for region in out_region[chrn]:
				f_out.write(chrn+'\t'+str(region[0])+'\t'+str(region[1])+'\n')


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Notice: this script is using for remove regions from choromosome with blast results")
		print("Usage: python "+sys.argv[0]+" <blast_results> <chr_len> <out_bed>")
		print("\t<blast_results> is a blast file list segmented with comma")
	else:
		blast_results = sys.argv[1].split(",")
		chr_len_file = sys.argv[2]
		out_bed = sys.argv[3]
		remove_region(blast_results, chr_len_file, out_bed)

