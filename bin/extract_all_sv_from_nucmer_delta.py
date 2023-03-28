#!/usr/bin/env python
import sys, os


def nucmer_extract(in_delta, out_pre):
	print("Running delta-filter")
	cmd = "delta-filter -gqr "+in_delta+" > "+in_delta+".filtered"
	print("Running command: "+cmd)
	os.system(cmd)
	print("Extracting")
	data_db = {}
	last_INDEL_pos = {}
	r_chr_len_db = {}
	q_chr_len_db = {}
	#sv_list = ['SNP', 'INDEL', 'JMP', 'INV', 'DUP', 'BRK']
	with os.popen("show-diff -H "+in_delta+".filtered", 'r') as f_in:
		for line in f_in:
			data = line.strip().split()
			if len(data) < 5:
				continue
			chrn = data[0]
			sv = data[1]
			if sv not in data_db:
				data_db[sv] = {}
			if chrn not in data_db[sv]:
				data_db[sv][chrn] = []
			sp = data[2]
			ep = data[3]
			data_db[sv][chrn].append([sp, ep])			
	
	data_db['SNP'] = {}
	data_db['INDEL'] = {}
	with os.popen("show-snps -ClrT "+in_delta+".filtered", 'r') as f_in:
		for line in f_in:
			data = line.strip().split()
			if len(data) == 0 or data[0].isdigit() == False:
				continue
			
			pos = int(data[0])
			r_chrn = data[-2]
			if r_chrn not in data_db['SNP']:
				data_db['SNP'][r_chrn] = []
			if r_chrn not in data_db['INDEL']:
				data_db['INDEL'][r_chrn] = []
			if data[1] != '.' and data[2] != '.':
				data_db['SNP'][r_chrn].append([data[0], data[1], data[2]])
			else:
				data_db['INDEL'][r_chrn].append([data[0], data[1], data[2]])
	
	print("Writing data")
	for sv in data_db:
		with open(out_pre+"."+sv+".txt", 'w') as fout:
			for chrn in sorted(data_db[sv]):
				for data in data_db[sv][chrn]:
					fout.write("%s\t%s\n"%(chrn, '\t'.join(data)))
	print("Success")


if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Usage: python "+sys.argv[0]+" <in_delta> <out_pre>")
	else:
		in_delta, out_pre = sys.argv[1:]
		nucmer_extract(in_delta, out_pre)
