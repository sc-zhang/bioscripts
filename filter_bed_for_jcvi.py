#!/usr/bin/env python
import sys


def filter_bed(in_bed, chrn, sp, ep, out_bed):
	sp = int(sp)
	ep = int(ep)
	with open(in_bed, 'r') as f_in:
		with open(out_bed, 'w') as f_out:
			for line in f_in:
				data = line.strip().split()
				cur_chrn = data[0]
				if cur_chrn != chrn:
					continue
				cur_sp = int(data[1])
				cur_ep = int(data[2])
				if cur_sp >= sp and cur_sp <= ep or cur_ep >= sp and cur_ep <= ep:
					f_out.write(line)


if __name__ == "__main__":
	if len(sys.argv) < 6:
		print("Usage: python "+sys.argv[0]+" <in_bed> <chr_name> <start_pos> <end_pos> <out_bed>")
	else:
		prog, in_bed, chrn, sp, ep, out_bed = sys.argv
		filter_bed(in_bed, chrn, sp, ep, out_bed)

