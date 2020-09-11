#!/usr/bin/env python
import sys, os


def generate_scripts(in_geno, in_list, out_dir):
	if os.path.isdir(out_dir) == False:
		os.mkdir(out_dir)
	
	if '/' not in in_geno:
		in_geno = os.path.join(os.getcwd(), in_geno)
	pop_db = []
	with open(in_list, 'r') as f_in:
		for line in f_in:
			data = line.strip().split()
			if data[0] == 'D':
				continue
			else:
				pop = [data[8], data[9], data[10], data[11]]
				if pop not in pop_db:
					pop_db.append(pop)
	
	for i in range(0, len(pop_db)):
		with open(out_dir+"/run_ABBA_"+str(i)+".sh", 'w') as f_out:
			f_out.write("#!/bin/bash\n#$ -cwd\n#$ -S /bin/bash\n#$ -j y\n#$ -pe mpi 10\n#$ -q all.q\n")
			f_out.write("python /public1/home/zsc/software/genomics_general/ABBABABAwindows.py -g "+in_geno+" -f phased -o output_"+"-".join(pop_db[i])+".csv -w 100000 -m 100 -s 100000 -P1 A -P2 B -P3 C -O D -T 10 --minData 0.5 --popsFile pops_"+"-".join(pop_db[i])+".txt --writeFailedWindows")

		with open(out_dir+"/pops_"+"-".join(pop_db[i])+".txt", 'w') as f_out:
			f_out.write(pop_db[i][0]+"\tA\n"+pop_db[i][1]+"\tB\n"+pop_db[i][2]+"\tC\n"+pop_db[i][3]+"\tD\n")


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Notice: script for generating scripts for ABBABABAwindows.py with ErrorCorr file")
		print("Usage: python "+sys.argv[0]+" <in_geno> <in_list> <out_dir>")
	else:
		in_geno = sys.argv[1]
		in_list = sys.argv[2]
		out_dir = sys.argv[3]
		generate_scripts(in_geno, in_list, out_dir)

