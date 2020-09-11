#!/usr/bin/python
import sys, os, shutil
import re


def make_group(db):
	new_db = {}
	for tig_name in db:
		max_name = ''
		max = 0
		for chrn in db[tig_name]:
			if max < db[tig_name][chrn]:
				max_name = chrn
		if chrn not in new_db:
			new_db[chrn] = []
		new_db[max_name].append(tig_name)
	return new_db


def split_group(blast_out, out_dir):
	if os.path.isdir(out_dir):
		shutil.rmtree(out_dir)
	os.mkdir(out_dir)
	
	print("Reading blast result")	
	blast_db = {}
	with open(blast_out, 'r') as f_in:
		for line in f_in:
			data = line.strip().split()
			if data[1][:3].lower() != 'chr':
				continue
			if data[1][3].isdigit() == False:
				continue
			tig_name = data[0].split("_")[0]
			chrn = data[1]
			if tig_name not in blast_db:
				blast_db[tig_name] = {}
			if chrn not in blast_db[tig_name]:
				blast_db[tig_name][chrn] = 0
			blast_db[tig_name][chrn] += 1
	print("Making groups")
	group_db = make_group(blast_db)
	
	f_db = {}
	print("Writing results")
	for chrn in group_db:
		f_db[chrn] = open(out_dir+"/"+chrn+".ordering", 'w')
		f_db[chrn].write("\n".join(group_db[chrn]))
		f_db[chrn].close()	
	
	print("Success")


if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Notice: group contigs with blast result")
		print("Usage python "+sys.argv[0]+" <blast_out> <out_dir>")
	else:
		blast_out = sys.argv[1]
		out_dir = sys.argv[2]
		split_group(blast_out, out_dir)

