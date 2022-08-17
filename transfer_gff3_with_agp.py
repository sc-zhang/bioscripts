#!/usr/bin/env python
import sys


def get_gene_ctg(chr_list, sp, ep):
	infos = []
	for info in chr_list:
		chsp = info[0]
		chep = info[1]
		ovlp = min(chep, ep)-max(chsp, sp)
		if ovlp>=0:
			infos.append([ovlp, info])
	if infos != []:
		return sorted(infos, reverse=True)[0][1]
	else:
		return []


def trans_anno(in_gff3, in_old_agp, in_new_agp, out_gff3):
	print("Reading gff3")
	gff3_db = {}
	with open(in_gff3, 'r') as fin:
		for line in fin:
			if line[0] == '#' or line.strip() == '':
				continue
			data = line.strip().split()
			if data[2] == 'gene':
				ID = data[8].split(';')[0].split('=')[-1]
				if ID not in gff3_db:
					gff3_db[ID] = []
			gff3_db[ID].append(line)

	print("Reading old agp")
	old_agp_db = {}
	with open(in_old_agp, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			if data[4] == 'U':
				continue
			chrn = data[0]
			sp = int(data[1])
			ep = int(data[2])
			direct = data[-1]
			tig = data[5]
			if chrn not in old_agp_db:
				old_agp_db[chrn] = []
			old_agp_db[chrn].append([sp, ep, tig, direct])
	
	for chrn in old_agp_db:
		old_agp_db[chrn] = sorted(old_agp_db[chrn])
	
	print("Reading new agp")
	new_agp_db = {}
	with open(in_new_agp, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			if data[4] == 'U':
				continue
			tig = data[5]
			if 'r' in tig:
				tig, idx = tig.split('r')
				if idx != '1':
					continue
			chrn = data[0]
			sp = int(data[1])
			ep = int(data[2])
			direct = data[-1]
			new_agp_db[tig] = [chrn, sp, ep, direct]
	
	print("Writing new gff3")
	with open(out_gff3, 'w') as fout:
		fout.write("###gff version 3\n")
		for id in sorted(gff3_db):
			for i in range(0, len(gff3_db[id])):
				data = gff3_db[id][i].split()
				if data[2] == 'gene':
					break
			chrn = data[0]
			#if chrn not in old_agp_db:
			#	continue
			sp = int(data[3])
			ep = int(data[4])
			match_ctg = get_gene_ctg(old_agp_db[chrn], sp, ep)
			if match_ctg == []:
				print(id, data)
			else:
				csp, cep, tig, tdir = match_ctg
				nchrn, nsp, nep, ndir = new_agp_db[tig]
				for line in gff3_db[id]:
					data = line.strip().split()
					gsp = int(data[3])
					gep = int(data[4])
					gdir = data[6]
					if tdir == '+':
						gts = gsp-csp+1
						gte = gep-csp+1
					else:
						gts = cep-gep+1
						gte = cep-gsp+1
					if gdir == tdir:
						gtd = '+'
					else:
						gtd = '-'
					if ndir == '+':
						gns = nsp+gts-1
						gne = nsp+gte-1
					else:
						gns = nep-gte+1
						gne = nep-gts+1
					if gtd == ndir:
						gnd = '+'
					else:
						gnd = '-'
					if gns <= 0 or gne <= 0:
						continue
					data[0] = nchrn
					data[3] = str(gns)
					data[4] = str(gne)
					data[6] = gnd
					fout.write('\t'.join(data)+'\n')
			fout.write('\n')
	print("Finished")


if __name__ == "__main__":
	if len(sys.argv) < 5:
		print("Usage: python "+sys.argv[0]+" <in_gff3> <in_old_agp> <in_new_agp> <out_gff3>")
	else:
		in_gff3, in_old_agp, in_new_agp, out_gff3 = sys.argv[1:]
		trans_anno(in_gff3, in_old_agp, in_new_agp, out_gff3)
