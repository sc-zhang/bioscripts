#!/usr/bin/python
import sys


def calc_gene_density(in_gene, in_bed, out_gene):
	region_db = {}
	win_length = 0
	with open(in_bed, 'r') as f_in:
		for line in f_in:
			if line[:3].lower() != 'chr':
				continue
			data = line.strip().split()
			if win_length == 0:
				win_length = int(data[2]) - int(data[1])
			if data[0] not in region_db:
				region_db[data[0]] = {}
				region_db[data[0]]['start'] = []
				region_db[data[0]]['count'] = []
			region_db[data[0]]['start'].append(int(data[1]))
			region_db[data[0]]['count'].append(0)
	with open(in_gene, 'r') as f_in:
		for line in f_in:
			if line[:3].lower() != 'chr':
				continue
			data = line.strip().split()
			chrn = data[0]
			s = int(data[1])
			for i in range(0, len(region_db[chrn]['start'])):
				if s < region_db[chrn]['start'][i]:
					break
			i -= 1
			region_db[chrn]['count'][i] += 1
	chr_list = sorted(region_db.keys())
	with open(out_gene, 'w') as f_out:
		for chrn in chr_list:
			for i in range(0, len(region_db[chrn]['start'])):
				f_out.write(chrn+'\t'+str(region_db[chrn]['start'][i])+'\t'+str(region_db[chrn]['start'][i]+win_length)+'\t'+str(region_db[chrn]['count'][i])+'\n')


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Notice: calc gene density counts with bed file")
		print("Usage: python "+sys.argv[0]+" <in_gene> <in_bed> <out_gene>")
	else:
		in_gene = sys.argv[1]
		in_bed = sys.argv[2]
		out_gene = sys.argv[3]
		calc_gene_density(in_gene, in_bed, out_gene)

