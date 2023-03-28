#!/usr/bin/env python
import sys


def merge_regions(regions):
	tmp_regions = []
	last_ep = 0
	for sp, ep in sorted(regions):
		if tmp_regions == []:
			tmp_regions.append(sp)
			last_ep = ep
		else:
			if sp > last_ep:
				tmp_regions.append(last_ep)
				tmp_regions.append(sp)
				last_ep = ep
			else:
				if ep > last_ep:
					last_ep = ep
	tmp_regions.append(last_ep)
	new_regions = []
	for i in range(0, len(tmp_regions)-1, 2):
		new_regions.append([tmp_regions[i], tmp_regions[i+1]])
	return new_regions


def calc_ovlp_ratio(regions, sp, ep):
	tmp_tes = []
	for rsp, rep in regions:
		ovlp = min(ep, rep)-max(sp, rsp)+1
		if ovlp <= 0:
			continue
		tmp_tes.append([max(sp, rsp), min(ep, rep)])
	
	ovlp_len = 0
	if len(tmp_tes) != 0:
		for msp, mep in merge_regions(tmp_tes):
			ovlp_len += mep-msp+1
	
	return ovlp_len*1.0/(ep-sp+1)


def calc_gene_ovlp_te(gene_gff3, TE_gffs, ovlp_stat):
	print("Loading TEs")
	TE_db = {}
	for te in TE_gffs.split(','):
		print("\tLoading: %s"%te)
		with open(te, 'r') as fin:
			for line in fin:
				if line[0] == '#':
					continue
				data = line.strip().split()
				tig = data[0]
				sp = int(data[3])
				ep = int(data[4])
				if sp > ep:
					sp, ep = ep, sp
				if tig not in TE_db:
					TE_db[tig] = []
				TE_db[tig].append([sp, ep])

	for tig in TE_db:
		TE_db[tig] = sorted(TE_db[tig])
	
	print("Reading gene gff3 and calculating overlaps")
	with open(gene_gff3, 'r') as fin:
		with open(ovlp_stat, 'w') as fout:
			for line in fin:
				if line[0] == '#':
					continue
				data = line.strip().split()
				if data[2] != 'gene':
					continue
				tig = data[0]
				gid = data[8].split(';')[0].split('=')[1]
				sp = int(data[3])
				ep = int(data[4])
				if sp > ep:
					sp, ep = ep, sp
				if tig not in TE_db:
					continue
				fout.write("%s\t%f\n"%(gid, 100*calc_ovlp_ratio(TE_db[tig], sp, ep)))
	print("Finished")


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python "+sys.argv[0]+" <gene_gff3> <TE_gffs> <ovlp_stat>")
	else:
		gene_gff3, TE_gffs, ovlp_stat = sys.argv[1:]
		calc_gene_ovlp_te(gene_gff3, TE_gffs, ovlp_stat)
