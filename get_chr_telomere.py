#!/usr/bin/env python
import sys
import copy


def get_chr_cen(in_ctg_cen, in_agp, bin_size, out_chr_cen):
	print("Loading agp")
	chr_len = {}
	ctg_on_chr = {}
	with open(in_agp, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			if data[4]!= 'W' or data[0] == data[5]:
				continue
			chrn = data[0]
			tig = data[5]
			sp = int(data[1])
			ep = int(data[2])
			direct = data[-1]
			if chrn not in chr_len:
				chr_len[chrn] = 0
			chr_len[chrn] = ep
			ctg_on_chr[tig] = [chrn, sp, ep, direct]
	
	telo_dist = {}
	for chrn in chr_len:
		telo_dist[chrn] = [ 0 for i in range(0, int(chr_len[chrn]/bin_size+2))]
	
	print("Loading telomere information with contig level")
	telo_db = {}
	with open(in_ctg_cen, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			if data[0] == 'seqid':
				continue
			tig = data[0]
			tsp = int(data[1])
			tep = int(data[2])
			if tig not in ctg_on_chr:
				continue
			chrn, csp, cep, direct = ctg_on_chr[tig]
			if direct == '+':
				nsp = csp+tsp-1
				nep = csp+tep-1
			else:
				nsp = cep-tep+1
				nep = cep-tsp+1
			if chrn not in telo_db:
				telo_db[chrn] = []
			telo_db[chrn].append([nsp, nep])
	
	for chrn in telo_db:
		telo_db[chrn] = sorted(telo_db[chrn])
		for sp, ep in telo_db[chrn]:
			sidx = int(round(sp/bin_size+0.51))
			eidx = int(round(sp/bin_size+0.51))
			for i in range(sidx, eidx+1):
				telo_dist[chrn][i] += 1
	
	print("Calculating telomere positions")
	res_db = {}
	for chrn in sorted(telo_dist):
		print("\tCalculating %s"%chrn)
		tmp_list = sorted(telo_dist[chrn], reverse=True)
		d = len(tmp_list)/100.0
		last_idx = -1
		telo_pre = []
		for i in range(0, len(telo_dist[chrn])):
			cur_val = telo_dist[chrn][i]
			if cur_val > 0:
				if last_idx == -1:
					telo_pre.append([i, i])
				else:
					if i <= last_idx+d:
						telo_pre[-1][1] = i
					else:
						telo_pre.append([i, i])
				last_idx = i
		max_pair = []
		sec_max_pair = []
		max_len = 0
		sec_max_len = 0
		for si, ei in telo_pre:
			if ei-si+1 > max_len:
				if max_len != 0:
					sec_max_len = max_len
					sec_max_pair = max_pair
				max_len = ei-si+1
				max_pair = [si, ei]
		if max_pair == []:
			ms = "NA"
			me = "NA"
			ml = "NA"
		else:
			ms = max_pair[0]*bin_size
			me = max_pair[1]*bin_size
			ml = me-ms+1
			ms = "%.2f"%(ms/1e6)
			me = "%.2f"%(me/1e6)
			ml = "%.2f"%(ml/1e6)
		if sec_max_pair == []:
			ss = "NA"
			se = "NA"
			sl = "NA"
		else:
			ss = sec_max_pair[0]*bin_size
			se = sec_max_pair[1]*bin_size
			sl = se-ss+1
			ss = "%.2f"%(ss/1e6)
			se = "%.2f"%(se/1e6)
			sl = "%.2f"%(sl/1e6)
		base_chr = chrn[:-1]
		hap = chrn[-1]
		if base_chr not in res_db:
			res_db[base_chr] = {}
		res_db[base_chr][hap] = [ms, me, ml, ss, se, sl]
	
	hap_list = sorted(res_db[base_chr])

	print("Writing result")
	with open(out_chr_cen, 'w') as fout:
		fout.write(",%s\n"%(',,,,'.join(hap_list)))
		for chrn in sorted(res_db):
			info = []
			for hap in sorted(res_db[chrn]):
				ms, me, ml, ss, se, sl = res_db[chrn][hap]
				info.append("%s-%s,%s,%s-%s,%s"%(ms, me, ml, ss, se, sl))
			fout.write("%s,%s\n"%(chrn, ','.join(info)))
	
	print("Finished")


if __name__ == "__main__":
	if len(sys.argv) < 5:
		print("Usage: python %s <in_ctg_telo> <in_agp> <bin_size> <out_chr_telo>"%sys.argv[0])
	else:
		in_ctg_cen, in_agp, bin_size, out_chr_cen = sys.argv[1:]
		if not bin_size[-1].isdigit():
			ratio = bin_size[-1].lower()
			base = float(bin_size[:-1])
			if ratio == 'g':
				ratio = 1e9
			elif ratio == 'm':
				ratio = 1e6
			else:
				ratio = 1e3
			bin_size = int(base*ratio)
		else:
			bin_size = int(bin_size)
		get_chr_cen(in_ctg_cen, in_agp, bin_size, out_chr_cen)
