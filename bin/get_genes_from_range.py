#!/usr/bin/python
import sys


def cmplist(a,b):
	dataa = a.strip().split('\t')
	datab = b.strip().split('\t')
	if long(dataa[0]) > long(datab[0]):
		return 1
	elif long(dataa[0]) < long(datab[0]):
		return -1
	else:
		return 0

def search(listGene, x, y, td):
	res = ""
	for line in listGene:
		data = line.split('\t')
		if (x >= long(data[0]) and x <= long(data[1])) or (y >= long(data[0]) and y <= long(data[1])) or (x <= long(data[0]) and y>= long(data[1])):
			if x >= long(data[0]):
				start_pos = x
			else:
				start_pos = long(data[0])
			if y <= long(data[1]):
				end_pos = y
			else:
				end_pos = long(data[1])
			if (end_pos-start_pos+1)*1.0/(y-x+1) >= td:
				res = res+data[2]+'\n'
	return res


def get_genes_from_range(f_gff, f_bed, f_out, td):
	gff = open(f_gff,'r')
	dictGene = {}
	for line in gff:
		data = line.strip().split('\t')
		if(len(data) > 3):
			if(data[2] == "gene"):
				if data[0] not in dictGene:
					dictGene[data[0]] = []
				dictGene[data[0]].append(data[3]+'\t'+data[4]+'\t'+data[8].split(';')[0][3:])
	gff.close()
	
	for key in dictGene:
		dictGene[key].sort(cmplist)
	bed = open(f_bed,'r')
	out = open(f_out,'w')
	for line in bed:
		data = line.strip().split('\t')
		if len(data) < 3 or line[0] == '#':
			continue
		ss = long(data[1])
		se = long(data[2])
		if data[0] in dictGene:
			res = search(dictGene[data[0]], ss, se, td)
			if res != "":
				rs = res.split('\n')
				out.write(line.strip())
				for r in rs:
					if r != "":
						out.write('\t'+r)
				out.write('\n')
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
