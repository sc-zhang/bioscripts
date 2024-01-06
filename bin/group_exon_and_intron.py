#!/usr/bin/env python
import sys
import time
import gzip


def search_pos(in_gff, in_vcf, out_file):
	pos_db = {}
	if in_gff.split('.')[-1] == 'gz':
		f_gff = gzip.open(in_gff, 'rt')
	else:
		f_gff = open(in_gff, 'r')

	for line in f_gff:
		if line[0] == '#' or line.strip() == '':
			continue
		data = line.strip().split()
		chrn = data[0]
		s_pos = int(data[3])
		e_pos = int(data[4])
		name = data[8].split(';')[0].split('=')[1]
		if "G" in name:
			name = '.'.join(name.split('.')[:2])
		else:
			name = name.split('.')[0]
		type = data[2]
		if chrn not in pos_db:
			pos_db[chrn] = {}
		if name not in pos_db[chrn]:
			pos_db[chrn][name] = {}
			pos_db[chrn][name]['gene'] = ()
			pos_db[chrn][name]['exon'] = []
		if type == 'gene':
			pos_db[chrn][name]['gene'] = (s_pos, e_pos)
		elif type == 'exon':
			pos_db[chrn][name]['exon'].append((s_pos, e_pos))
		else:
			continue
	f_gff.close()
	
	if in_vcf.split('.') == 'gz':
		f_vcf = gzip.open(in_vcf, 'rt')
	else:
		f_vcf = open(in_vcf, 'r')
	
	if out_file.split('.') == 'gz':
		f_out = gzip.open(out_file, 'wt')
	else:
		f_out = open(out_file, 'w')

	for line in f_vcf:
		if line[0] == '#' or line.strip() == '':
			continue
		data = line.strip().split()
		chrn = data[0]
		pos = int(data[1])
		if chrn not in pos_db:
			continue
		is_found = False
		for name in pos_db[chrn]:
			s_pos, e_pos = pos_db[chrn][name]['gene']
			if s_pos <= pos and pos <= e_pos:
				is_found = True
				break
		if is_found:
			is_exon = False
			for (s_pos, e_pos) in pos_db[chrn][name]['exon']:
				if s_pos <= pos and pos <= e_pos:
					is_exon = True
					break
			if is_exon:
				f_out.write(chrn+"\t"+name+"\t"+str(pos)+"\t"+"exon\n")
			else:
				f_out.write(chrn+"\t"+name+"\t"+str(pos)+"\t"+"intron\n")
	f_vcf.close()
	f_out.close()


if __name__ == "__main__":
	if len(sys.argv)<4:
		print("Notice: script use for determining pos to exon or intro in vcf file base on gff file")
		print("Usage python "+sys.argv[0]+" <input_gff> <input_vcf> <output_file>")
	else:
		s_time = time.time()
		in_gff = sys.argv[1]
		in_vcf = sys.argv[2]
		out_file = sys.argv[3]
		search_pos(in_gff, in_vcf, out_file)
		e_time = time.time()
		print("cost time " + str(e_time-s_time))

