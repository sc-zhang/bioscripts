#!/usr/bin/python
import sys, os, shutil


def parse_gff(gff_file):
	gff_db = {}
	with open(gff_file, "r") as f_gff:
		for lines in f_gff:
			if lines[0] == "#":
				continue
			data = lines.strip().split()
			if data[2] != "gene":
				continue
			g_pos = data[0]
			g_name = data[8].split(";")[1].split("=")[1]
			if g_name not in gff_db:
				gff_db[g_name] = {}
			if g_pos not in gff_db[g_name]:
				gff_db[g_name][g_pos] = 0
			gff_db[g_name][g_pos] += 1
	return gff_db


def max_map(gff_db):
	for g_name in gff_db:
		max_pos = ""
		for g_pos in gff_db[g_name]:
			if max_pos == "":
				max_pos = g_pos
			else:
				if gff_db[g_name][g_pos] > gff_db[g_name][max_pos]:
					max_pos = g_pos
		gff_db[g_name] = max_pos
	return gff_db

def merge_map(ref_db, query_db, map_db):
	query_map = {}
	for g_name in ref_db:
		if g_name in map_db:
			if map_db[g_name] in query_db:
				q_name = query_db[map_db[g_name]]
				if q_name not in query_map:
					query_map[q_name] = {}
				if ref_db[g_name] not in query_map[q_name]:
					query_map[q_name][ref_db[g_name]] = 0
				query_map[q_name][ref_db[g_name]] += 1
	query_map = max_map(query_map)
	merge_db = {}
	for g_query in query_map:
		if query_map[g_query] not in merge_db:
			merge_db[query_map[g_query]] = []
		merge_db[query_map[g_query]].append(g_query)
	return merge_db;



def split_group(blast_result, ref_gff_file, query_gff_file, out_dir):
	if os.path.isdir(out_dir):
		shutil.rmtree(out_dir)
	os.mkdir(out_dir)
	
	print("Read blast result")
	blast_db = {}
	with open(blast_result, "r") as f_blast:
		for lines in f_blast:
			data = lines.strip().split()
			blast_db[data[1]] = data[0]
	
	print("Read gff file")
	
	ref_chr_db = parse_gff(ref_gff_file)
	query_ctg_db = parse_gff(query_gff_file)
		
	print("Mapping")
	ref_chr_db = max_map(ref_chr_db)
	query_ctg_db = max_map(query_ctg_db)
	
	print("Merging")
	merge_db = merge_map(ref_chr_db, query_ctg_db, blast_db)
	
	print("Writing")
	for chrn in merge_db:
		if chrn[:3].lower() != "chr" or chrn[3].isdigit() == False:
			continue
		f_out_name = out_dir + "/group" + chrn[-1] + ".ordering"
		with open(f_out_name, "w") as f_out:
			for ctg in merge_db[chrn]:
				if ctg[:3] != "tig":
					continue
				f_out.write(ctg+"\n")
				
	print("Success")


if __name__ == "__main__":
	if len(sys.argv) > 4:
		blast_result = sys.argv[1]
		ref_gff_file = sys.argv[2]
		query_gff_file = sys.argv[3]
		out_dir = sys.argv[4]
		split_group(blast_result, ref_gff_file, query_gff_file, out_dir)
	else:
		print("Notice: script for spliting query contigs into different groups of chromosomes with blast result of contig files")
		print("Usage: python "+sys.argv[0]+" <blast_result> <ref_gff> <query_gff> <out_dir>")

