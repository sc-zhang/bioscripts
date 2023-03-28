#!/usr/bin/env python
import sys


def convert_gbff_to_fasta(in_gbff, out_fa):
	print("Converting")
	with open(in_gbff, 'r') as fin:
		with open(out_fa, 'w') as fout:
			cnt = 0
			err_cnt = 0
			for line in fin:
				data = line.strip().split()
				if data[0] == 'LOCUS':
					cnt += 1
					gn = data[1]
					gn_len = int(data[2])
					seq_len = 0
					fout.write(">%s\n"%gn)
					is_write = False
				elif data[0] == 'ORIGIN':
					is_write = True
				elif data[0] == '//':
					if gn_len != seq_len:
						err_cnt += 1
						print("\tERROR: %s Comment length: %sbp, current length: %dbp"%(gn, gn_len, seq_len))
					is_write = False
				else:
					if is_write:
						seq = ''.join(data[1:])
						seq_len += len(seq)
						fout.write(seq+'\n')
	print("Total convert %d, error count %d"%(cnt, err_cnt))
	print("Finished")
	

if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Usage: python %s <in_gbff> <out_fasta>"%sys.argv[0])
	else:
		in_gbff, out_fa = sys.argv[1:]
		convert_gbff_to_fasta(in_gbff, out_fa)
