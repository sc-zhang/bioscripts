#!/usr/bin/env python
import sys, os


def nucmer_extract(r_fa, q_fa, out_pre, ts):
	print("Running nucmer")
	cmd = "nucmer -p "+out_pre+" "+r_fa+" "+q_fa+" -t "+ts
	print("Running command: "+cmd)
	os.system(cmd)

	print("Running delta-filter")
	cmd = "delta-filter -gqr "+out_pre+".delta > "+out_pre+".filtered"
	print("Running command: "+cmd)
	os.system(cmd)

	print("Extracting")
	data_db = {}
	last_INDEL_pos = {}
	r_chr_len_db = {}
	q_chr_len_db = {}
	#sv_list = ['SNP', 'INDEL', 'JMP', 'INV', 'DUP', 'BRK']
	with os.popen("show-diff -H "+out_pre+".filtered", 'r') as f_in:
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
	with os.popen("show-snps -ClrT "+out_pre+".filtered", 'r') as f_in:
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
	with open(out_pre+".sv.txt", 'w') as fall:
		for sv in data_db:
			with open(out_pre+"."+sv+".txt", 'w') as fout:
				for chrn in sorted(data_db[sv]):
					for data in data_db[sv][chrn]:
						if sv != "SNP" and sv != "INDEL":
							fall.write("%s\t%s\t%s\n"%(chrn, '\t'.join(data), sv))
						fout.write("%s\t%s\n"%(chrn, '\t'.join(data)))
	print("Success")


if __name__ == "__main__":
	if len(sys.argv) < 5:
		print("Usage: python "+sys.argv[0]+" <ref_fasta> <query_fasta> <out_pre> <threads>")
	else:
		r_fa, q_fa, out_pre, ts = sys.argv[1:]
		nucmer_extract(r_fa, q_fa, out_pre, ts)
