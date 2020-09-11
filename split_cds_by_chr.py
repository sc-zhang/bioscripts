#!/usr/bin/python
import sys


def split_cds(bed_list, cds_file):
	cds_db = {}
	with open(cds_file, 'r') as f_in:
		id = ''
		seq = ''
		for line in f_in:
			line = line.strip()
			if line[0] == '>':
				if seq != '':
					cds_db[id] = seq
				id = line
				seq = ''
			else:
				seq += line
		cds_db[id] = seq
	
	read_files = {}
	write_files = {}
	with open(bed_list, 'r') as f_in:
		for line in f_in:
			line = line.strip()
			if line != '' and line not in read_files:
				read_files[line] = open(line, 'r')
				write_files[line] = open(line.replace("bed", "cds"), 'w')

	for fn in read_files:
		for line in read_files[fn]:
			data = line.strip().split()
			id = ">"+data[3]
			if id not in cds_db:
				continue
			write_files[fn].write(id+"\n"+cds_db[id]+"\n")
		read_files[fn].close()
		write_files[fn].close()


if __name__ == "__main__":
	if len(sys.argv) < 2:
		print("Usage: python "+sys.argv[0]+" <bed_list> <cds_file>")
	else:
		bed_list = sys.argv[1]
		cds_file = sys.argv[2]
		split_cds(bed_list, cds_file)

