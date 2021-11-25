#!/usr/bin/env python
import sys


def generate_new_id_by_gff(chr_pre, in_gff):
	id_db = {}
	gff_db = {}
	is_first_id = True
	with open(in_gff, 'r') as f_in:
		for line in f_in:
			data = line.strip().split()
			if line[0] == '#' or len(data) < 9:
				continue
			if data[2] == 'gene':
				chrn = data[0]
			#if chrn[:3].lower() != 'chr':
			#	continue
				s_p = int(data[3])
				id = data[8].split(";")[0].split("=")[1]
				if is_first_id:
					print("Check ID: %s"%id)
					is_first_id = False
				if chrn not in id_db:
					id_db[chrn] = []
				id_db[chrn].append([s_p, id])
				gff_db[id] = []
			gff_db[id].append(line)
				
	new_id_db = {}
	chr_base = {}
	ordered_id = []
	with open("rename_list.txt", 'w') as f_out:
		tig_base = 10
		for chrn in sorted(id_db):
			base = 10
			for info in sorted(id_db[chrn]):
				if chrn[:3].lower() == 'chr':
					idx = int(chrn[3:])
					new_id = chr_pre+".%02dG%07d"%(idx, base)
					base += 10
				else:
					new_id = chr_pre+".%08d"%(tig_base)
					tig_base += 10
				#	new_id = info[1].replace('.', '').replace('G', '')
				f_out.write("%s\t%d\t%s\t%s\n"%(chrn, info[0], info[1], new_id))
				ordered_id.append(info[1])
				new_id_db[info[1]] = new_id
	return new_id_db, gff_db, ordered_id
			

def rename_id(chr_pre, in_gff, out_gff, in_fastas, out_fastas):
	print("Generating rename list")
	rename_id_db, gff_db, ordered_id = generate_new_id_by_gff(chr_pre, in_gff)

	print("Dealing gff")
	with open(out_gff, 'w') as fout:
		fout.write("###gff version 3\n")
		for id in ordered_id:
			for line in gff_db[id]:
				fout.write(line.replace(id, rename_id_db[id]))
			fout.write("\n")				

	print("Dealing fasta")
	in_fasta_list = in_fastas.split(',')
	out_fasta_list = out_fastas.split(',')
	for i in range(0, len(in_fasta_list)):
		print("\tDealing %s"%in_fasta_list[i])
		with open(in_fasta_list[i], 'r') as f_in:
			with open(out_fasta_list[i], 'w') as f_out:
				suf = in_fasta_list[i].split('.')[-1].lower()
				if suf == "fasta" or suf == "fa":
					for line in f_in:
						if line[0] == '>':
							id = line.strip()[1:]
							if id in rename_id_db:
								line = line.replace(id, rename_id_db[id])
						f_out.write(line)
				else:
					for line in f_in:
						id = line.strip().split()[0]
						if id in rename_id_db:
							f_out.write(line.replace(id, rename_id_db[id]))
						else:
							f_out.write(line)
	print("Finished")


if __name__ == "__main__":
	if len(sys.argv) < 5:
		print("Usage: python "+sys.argv[0]+" <chr_prefix> <in_gff> <out_gff> <in_fasta> <out_fasta>")
		print("Notice: sort and rename id with in_gff, and rename them in fasta files")
		print("Example: python "+sys.argv[0]+" CB5 in.gff out.gff 1.fasta,2.fasta 1.new.fasta,2.new.fasta")
	else:
		chr_pre, in_gff, out_gff, in_fastas, out_fastas = sys.argv[1:]
		rename_id(chr_pre, in_gff, out_gff, in_fastas, out_fastas)
