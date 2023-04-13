#!/usr/bin/env python
import sys
import gzip


def subVCF(in_vcf, in_list, out_vcf, missing_rate):
	sp_list = []
	
	with open(in_list, 'r') as f_in:
		for line in f_in:
			if line.strip() != '':
				sp_list.append(line.strip())
	
	if in_vcf.split('.')[-1] == 'gz':
		f_in = gzip.open(in_vcf, 'rt')
	else:
		f_in = open(in_vcf, 'r')
	

	if out_vcf.split('.')[-1] == 'gz':
		f_out = gzip.open(out_vcf, 'wt')
	else:
		f_out = open(out_vcf, 'w')
	
	col_db = []
	for line in f_in:
		data = line.strip().split()
		pub_info = data[0:9]
		if data[0][0] == "#":
			for i in range(9, len(data)):
				if data[i] in sp_list:
					col_db.append(i)
		
		cnt_mis = 0
		out_str = ''
		for i in col_db:
			out_str += "\t"+data[i]
			if data[i] == './.' or data[i] == '.|.':
				cnt_mis += 1
		if len(col_db) > 0 and cnt_mis*1.0/len(col_db) > missing_rate:
			continue
		f_out.write("\t".join(pub_info))
		f_out.write(out_str+"\n")

	f_in.close()
	f_out.close()


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Notice: script to extract vcf file with list file, default missing rate 0.4")
		print("Usage: python "+sys.argv[0]+" <in_vcf> <in_list> <out_vcf> [<missing_rate>]")
	else:
		in_vcf = sys.argv[1]
		in_list = sys.argv[2]
		out_vcf = sys.argv[3]
		if len(sys.argv) == 5:
			missing_rate = float(sys.argv[4])
		else:
			missing_rate = 0.4
		subVCF(in_vcf, in_list, out_vcf, missing_rate)

