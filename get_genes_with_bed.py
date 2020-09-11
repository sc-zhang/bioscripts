#!/usr/bin/python
import sys
import re


def search(listGene, x, y, td):
	res = []
	for data in listGene:
		spos_r, epos_r, gn, cn = data
		ovlp = min(y, epos_r)-max(x, spos_r)
		if ovlp*1.0/(y-x+1) >= td:
			res.append(gn)
	return res


def get_genes_from_range(f_gff, f_bed, f_out, td):
	gff = open(f_gff,'r')
	dictGene = {}
	for line in gff:
		if line[0] == '#' or line.strip() == '':
			continue
		data = line.strip().split()
		if data[2] != "gene":
			continue
		if data[0] not in dictGene:
			dictGene[data[0]] = []
		gene_name = re.findall(r'Name=(.*?)$', data[8])[0]
		dictGene[data[0]].append([int(data[3]), int(data[4]), gene_name, data[0]])
	gff.close()
	
	
	out_infos = []
	for key in dictGene:
		dictGene[key] = sorted(dictGene[key])
	bed = open(f_bed,'r')
	out = open(f_out,'w')
	for line in bed:
		data = line.strip().split('\t')
		ss = int(data[1])
		se = int(data[2])
		if ss > se:
			ss, se = se, ss
		if data[0] in dictGene:
			res = search(dictGene[data[0]], ss, se, td)
			if res != []:
				out.write(line.strip()+'\t'+'\t'.join(res)+"\n")
	out.close()
	bed.close()


if __name__ == "__main__":
	if len(sys.argv) < 5:
		print("Notice: script for get genes in range from bed file")
		print("Usage: python "+sys.argv[0]+" <gff3_file> <bed_file> <output_file> <threshold>")
	else:
		f_gff = sys.argv[1]
		f_bed = sys.argv[2]
		f_out = sys.argv[3]
		td = float(sys.argv[4])/100.0
		get_genes_from_range(f_gff, f_bed, f_out, td)
