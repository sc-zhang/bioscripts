#!/usr/bin/env python
import sys, os
import gzip
import multiprocessing
import time
import gc


def read_fastq(in_fq):
	fq_db = {}
	id_list = []
	if in_fq[-3:].lower() == '.gz':
		f_fq = gzip.open(in_fq, 'rt')
	else:
		f_fq = open(in_fq, 'r')
	i = 0
	for line in f_fq:
		if i%4 == 0:
			id = line[1:]
			id_list.append(id)
			fq_db[id] = []
			i = 0
		i += 1
		fq_db[id].append(line.strip())
	f_fq.close()
	return id_list, fq_db


def generate_fragments(id_list, fq_db, out_fq, t_n):
	if out_fq[-3:].lower() == '.gz':
		f_out = gzip.open(out_fq+'_'+str(t_n), 'wt')
	else:
		f_out = open(out_fq+'_'+str(t_n), 'w')
	for id in id_list:
		read_length = len(fq_db[id][1])
		frag_cnt = int(read_length/36)
		sp = int((read_length%36)*2.0/3.0)
		for i in range(0, frag_cnt):
			if i*36+sp >= read_length or (i+1)*36+sp > read_length:
				continue
			f_out.write(fq_db[id][0]+'_'+str(i)+'\n')
			f_out.write(fq_db[id][1][i*36+sp: (i+1)*36+sp]+'\n')
			if fq_db[id][2] == '+':
				f_out.write('+\n')
			else:
				f_out.write(fq_db[id][2]+'_'+str(i)+'\n')
			f_out.write(fq_db[id][3][i*36+sp: (i+1)*36+sp]+'\n')
			#f_out.flush()
	f_out.close()


def fragment_reads(in_fq, out_fq, ts):
	start_t = time.time()
	print("Reading Fastq")
	id_list, fq_db = read_fastq(in_fq)
	task_per_process = int(len(id_list)/ts)
	task_list = []
	print("Generating fragments")
	
	for i in range(0, ts):
		if i < ts-1:
			t = multiprocessing.Process(target=generate_fragments, args=(id_list[i*task_per_process: (i+1)*task_per_process], fq_db, out_fq, i))
		else:
			t = multiprocessing.Process(target=generate_fragments, args=(id_list[i*task_per_process:], fq_db, out_fq, i))
		task_list.append(t)
	for t in task_list:
		t.start()
	for t in task_list:
		t.join()
	
	del id_list, fq_db
	gc.collect()
	print("Merging")
	if os.path.exists(out_fq):
		os.remove(out_fq)
	for i in range(0, ts):
		os.system("cat "+out_fq+"_"+str(i)+" >> "+out_fq)
		os.remove(out_fq+"_"+str(i))
	print("Success")
	end_t = time.time()
	print("Total cost time %.2fs"%(end_t-start_t))
	

if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python "+sys.argv[0]+" <in_fastq> <out_fastq> <threads>")
	else:
		prog, in_fq, out_fq, ts = sys.argv
		ts = int(ts)
		fragment_reads(in_fq, out_fq, ts)
