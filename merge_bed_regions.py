#!/usr/bin/env python
import sys

def merge_reions(region_db, md):
	new_regions = {}
	for chrn in region_db:
		new_regions[chrn] = []
		tmp_list = []
		tmp_gns = [[]]
		for region in sorted(region_db[chrn]):
			s = region[0]
			e = region[1]
			gn = region[2]
			if len(tmp_list) == 0:
				tmp_list.append(s)
				laste = e
			else:
				if s > laste+md:
					tmp_list.append(laste)
					tmp_list.append(s)
					tmp_gns.append([])
					laste = e
				else:
					if e > laste:
						laste = e
			tmp_gns[-1].extend(gn)
		tmp_list.append(laste)
		for i in range(0, len(tmp_list), 2):
			s = tmp_list[i]
			e = tmp_list[i+1]
			gns = ''.join(tmp_gns[int(i/2)])
			gns = list(set(gns.split(',')))
			new_gns = []
			for gn in gns:
				if gn != '':
					new_gns.append(gn)
			gns = ','.join(new_gns)
			new_regions[chrn].append([s, e, gns])
	return new_regions


def read_bed(in_bed):
	bed_db = {}
	with open(in_bed, 'r') as f_in:
		for line in f_in:
			if line.strip() == '':
				continue
			data = line.strip().split()
			chrn = data[0]
			sp = int(data[1])
			ep = int(data[2])
			if len(data) > 3:
				gn = data[-1]
			else:
				gn = ''
			if sp > ep:
				temp = sp
				sp = ep
				ep = temp
			if chrn not in bed_db:
				bed_db[chrn] = []
			bed_db[chrn].append([sp, ep, gn])
	return bed_db


def merge_regions_in_bed(in_bed, out_bed, md):
	ori_regions = read_bed(in_bed)
	new_regions = merge_reions(ori_regions, md)
	with open(out_bed, 'w') as f_out:
		for chrn in sorted(new_regions):
			for region in new_regions[chrn]:
				f_out.write("%s\t%d\t%d\t%s\n"%(chrn, region[0], region[1], region[2]))


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python "+sys.argv[0]+" <in_bed> <out_bed> <max_distance>")
	else:
		in_bed, out_bed, md = sys.argv[1:]
		md = int(md)
		merge_regions_in_bed(in_bed, out_bed, md)
