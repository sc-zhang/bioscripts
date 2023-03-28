#!/usr/bin/env python
import sys
import gzip
import time


def quick_extract_reads(in_fx, in_li, out_fx):
	print("\033[32m%s\033[0m Starting"%(time.strftime('[%H:%M:%S]',time.localtime(time.time()))))
	read_db = {}
	with open(in_li, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			read_db[data[0]] = ''
	
	fn = in_fx.split('.')
	if fn[-1].lower() == 'gz':
		fin = gzip.open(in_fx, 'rt')
	else:
		fin = open(in_fx, 'r')
	
	fn = out_fx.split('.')
	if fn[-1].lower() == 'gz':
		fout = gzip.open(out_fx, 'wt')
	else:
		fout = open(out_fx, 'w')
		
	is_write = False
	cnt = 0
	for line in fin:
		if cnt%2==0 and (line[0] == '>' or line[0] == '@'):
			id = line.strip().split()[0][1:]
			if id in read_db:
				is_write = True
				fout.write(line)
			else:
				is_write = False
		else:
			if is_write:
				fout.write(line)
		cnt += 1
	fin.close()
	fout.close()
	print("\033[32m%s\033[0m Finished"%(time.strftime('[%H:%M:%S]',time.localtime(time.time()))))


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python "+sys.argv[0]+" <in_fastx|gz> <in_list> <out_fastx|gz>")
	else:
		in_fx, in_li, out_fx = sys.argv[1:]
		quick_extract_reads(in_fx, in_li, out_fx)

