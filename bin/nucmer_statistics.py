#!/usr/bin/env python
import sys, os


def nucmer_statistics(r_fa, q_fa, out_pre, ts):
	print("Running nucmer")
	cmd = "nucmer -p "+out_pre+" "+r_fa+" "+q_fa+" -t "+ts
	print("Running command: "+cmd)
	os.system(cmd)
	
	print("Running delta-filter")
	cmd = "delta-filter -gqr "+out_pre+".delta > "+out_pre+".filtered"
	print("Running command: "+cmd)
	os.system(cmd)
	print("Statisitcs")
	data_db = {}
	last_INDEL_pos = {}
	r_chr_len_db = {}
	q_chr_len_db = {}
	sv_list = ['SNP', 'INDEL']
	with os.popen("show-snps -ClrT "+out_pre+".filtered", 'r') as f_in:
		for line in f_in:
			data = line.strip().split()
			if len(data) == 0 or data[0].isdigit() == False:
				continue
			
			pos = int(data[0])
			r_chrn = data[-2]
			q_chrn = data[-1]
			
			if r_chrn not in data_db:
				data_db[r_chrn] = {'SNP': {'count': 0, 'size': 0}, 'INDEL': {'count': 0, 'size': 0}}
				last_INDEL_pos[r_chrn] = 0
				r_chr_len_db[r_chrn] = int(data[6])
			if q_chrn not in q_chr_len_db:
				q_chr_len_db[q_chrn] = int(data[7])
			
			if data[1] != '.' and data[2] != '.':
				data_db[r_chrn]['SNP']['count'] += 1
				data_db[r_chrn]['SNP']['size'] += 1
				last_INDEL_pos[r_chrn] = 0
			else:
				if pos - last_INDEL_pos[r_chrn] > 1:
					data_db[r_chrn]['INDEL']['count'] += 1
				data_db[r_chrn]['INDEL']['size'] += 1
				last_INDEL_pos[r_chrn] = pos
	
	total_size = {}
	total_count = {}
	for chrn in data_db:
		total_size[chrn] = 0
		total_count[chrn] = 0
		for type in data_db[chrn]:
			total_size[chrn] += data_db[chrn][type]['size']
			total_count[chrn] += data_db[chrn][type]['count']
	
	with open(out_pre+".statistics", 'w') as f_out:
		f_out.write("Type\tNumber\tSize\tSize/TotalVar\tSize/ChrSize\n")
		for chrn in data_db:
			for type in sv_list:
				if type not in data_db[chrn]:
					continue
				f_out.write("%s\t%d\t%d\t%.4f\t%.4f\n"%(type, data_db[chrn][type]['count'], data_db[chrn][type]['size'], data_db[chrn][type]['size']*1.0/total_size[chrn], data_db[chrn][type]['size']*2.0/(r_chr_len_db[chrn]+q_chr_len_db[chrn])))
			f_out.write("%s\t%d\t%d\t%.4f\t%.4f\n"%("Total", total_count[chrn], total_size[chrn], total_size[chrn]*1.0/total_size[chrn], total_size[chrn]*2.0/(r_chr_len_db[chrn]+q_chr_len_db[chrn])))
	print("Success")


if __name__ == "__main__":
	if len(sys.argv) < 5:
		print("Usage: python "+sys.argv[0]+" <ref_fasta> <query_fasta> <out_pre> <threads>")
	else:
		r_fa, q_fa, out_pre, ts = sys.argv[1:]
		nucmer_statistics(r_fa, q_fa, out_pre, ts)
