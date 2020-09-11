#!/usr/bin/env python
import sys, os


def filter_blast(blast_file, out_file, t_i, t_m):
	blast_db = {}
	with open(blast_file, 'r') as f_in:
		for line in f_in:
			data = line.strip().split()
			chrn = data[1]
			iden = float(data[2])
			s_pos_1 = int(data[6])
			e_pos_1 = int(data[7])
			s_pos_2 = int(data[8])
			e_pos_2 = int(data[9])
			if iden < t_i:
				continue
			if s_pos_1 == s_pos_2 and e_pos_1 == e_pos_2:
				continue
			if s_pos_2 > e_pos_2:
				tmp = s_pos_2
				s_pos_2 = e_pos_2
				e_pos_2 = tmp
			if e_pos_2 - s_pos_2 < t_m:
				continue
			if chrn not in blast_db:
				blast_db[chrn] = []
			if [s_pos_2, e_pos_2] not in blast_db[chrn]:
				blast_db[chrn].append([s_pos_2, e_pos_2])
	with open(out_file, 'w') as f_out:
		for chrn in sorted(blast_db.keys()):
			for region in sorted(blast_db[chrn]):
				f_out.write(chrn+'\t'+str(region[0])+'\t'+str(region[1])+'\n')


def reshape(in_data, out_data):
	data_db = {}
	max_length = 0
	max_chr = ''
	with open(in_data, 'r') as f_in:
		for line in f_in:
			data = line.strip().split()
			chrn = data[0]
			pos = int(data[1])
			value = data[-1]
			if chrn[:3].lower() != 'chr':
				continue
			if len(chrn) == 4:
				chrn = chrn[:3]+'0'+chrn[-1]
			if chrn not in data_db:
				data_db[chrn] = []
			data_db[chrn].append([pos, value])

	for chrn in data_db:
		if len(data_db[chrn]) > max_length:
			max_length = len(data_db[chrn])
			max_chr = chrn
	
	for chrn in data_db:
		curr_len = len(data_db[chrn])
		if curr_len < max_length:
			for i in range(curr_len, max_length):
				data_db[chrn].append([data_db[max_chr][i][0], 'nan'])
	
	new_data = {}
	for chrn in data_db:
		for value in data_db[chrn]:
			if value[0] not in new_data:
				new_data[value[0]] = {}
			new_data[value[0]][chrn] = value[1]
	
	head = []
	for pos in new_data:
		for chrn in sorted(new_data[pos].keys()):
			head.append(chrn)
		break

	with open(out_data, 'w') as f_out:
		f_out.write('\t'+'\t'.join(head)+'\n')
		for pos in sorted(new_data.keys()):
			f_out.write(str(pos))
			for chrn in sorted(new_data[pos].keys()):
				f_out.write('\t'+new_data[pos][chrn])
			f_out.write('\n')


def draw_heatmap_R(in_data, out_name, rs):
	current_path = os.getcwd()
	script = os.path.join(current_path, out_name+'_draw.R')
	with open(script, 'w') as f_out:
		f_out.write("setwd(\""+current_path+"\")\n")
		f_out.write("data<-read.table(\""+in_data+"\", header = TRUE)\n")
		f_out.write("cr<-c(0:length(colnames(data))-1)\n")
		f_out.write("library(\"pheatmap\")\n")
		f_out.write("pheatmap(data, cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = FALSE, gaps_col = cr, filename = \""+out_name+".pdf\")\n")
	os.system("Rscript "+out_name+"_draw.R")

	script = os.path.join(current_path, out_name+'_draw_with_label.R')
	with open(script, 'w') as f_out:
		f_out.write("setwd(\""+current_path+"\")\n")
		f_out.write("data<-read.table(\""+in_data+"\", header = TRUE)\n")
		f_out.write("cr<-c(0:length(colnames(data))-1)\n")
		f_out.write("row_name<-array()\n")
		f_out.write("for(i in 1:length(rownames(data))){if(i%%==1){row_name[i]<-rownames(data)[i]}else{row_name[i]<-\"\"}}\n")
		f_out.write("library(\"pheatmap\")\n")
		f_out.write("pheatmap(data, cluster_cols = FALSE, cluster_rows = FALSE, labels_row = row_name, gaps_col = cr, filename = \""+out_name+"_label.pdf\")\n")
	os.system("Rscript "+out_name+"_draw_with_label.R")


def blast2heatmap(ref_fasta, blast_file, ws, out_name, t_i, t_m):
	print("Filter blast result")
	filter_blast(blast_file, "01_"+out_name+'_filter.bed', t_i, t_m)
	print("Generate windows")
	os.system("sliding_window.pl -i "+ref_fasta+" -o "+"02_"+out_name+"_win.bed -w "+ws+" -s "+ws)
	print("Coverage")
	os.system("bedtools coverage -a "+"02_"+out_name+"_win.bed -b "+"01_"+out_name+"_filter.bed > "+"03_"+out_name+"_cover.txt")
	print("Reshape data")
	reshape("03_"+out_name+"_cover.txt", "04_"+out_name+"_result.txt")
	print("Draw heatmap")
	draw_heatmap_R("04_"+out_name+"_result.txt", "05_"+out_name)


if __name__ == "__main__":
	if len(sys.argv) < 7:
		print("Usage: python "+sys.argv[0]+" <ref_fasta> <blast_file> <window_size> <out_name> <threshold_identify> <threshold_match>")
	else:
		proc, ref_fasta, blast_file, ws, out_name, t_i, t_m = sys.argv
		blast2heatmap(ref_fasta, blast_file, ws, out_name, float(t_i), float(t_m))
