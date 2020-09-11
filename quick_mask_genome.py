#!/usr/bin/env python
import sys, os
import multiprocessing


def read_bed(in_bed, th):
	win_regions = {}
	with open(in_bed, 'r') as f_in:
		for line in f_in:
			if line.strip() == '':
				continue
			data = line.strip().split()
			if int(data[3]) > 150:
				chrn = data[0]
				sp = int(data[1])-1
				ep = int(data[2])-1
				if chrn not in win_regions:
					win_regions[chrn] = []
				win_regions[chrn].append([sp, ep])
	return win_regions


def merge_regions(ori_regions):
	new_regions = {}
	for chrn in ori_regions:
		new_regions[chrn] = []
		tmp_regions = []
		last_e = 0
		for region in sorted(ori_regions[chrn]):
			sr = region[0]
			er = region[1]
			if last_e == 0:
				tmp_regions.append(sr)
				last_e = er
			if sr > last_e:
				tmp_regions.append(last_e)
				tmp_regions.append(sr)
				last_e = er
			else:
				if er > last_e:
					last_e = er
		tmp_regions.append(last_e)
		for i in range(0, len(tmp_regions), 2):
			new_regions[chrn].append([tmp_regions[i], tmp_regions[i+1]])
	return new_regions


def read_fasta(in_fa):
	seq_db = {}
	seq_id_list = []
	with open(in_fa, 'r') as f_in:
		id = ''
		seq = ''
		for line in f_in:
			if line[0] == '>':
				if seq != '':
					seq_db[id] = seq
				id = line.strip()[1:]
				seq_id_list.append(id)
				seq = ''
			else:
				seq += line.strip()
		seq_db[id] = seq
	return seq_db, seq_id_list


def mask_fasta(id_list, seq_db, win_regions):
	for id in id_list:
		new_seq = ''
		sp = 0
		if id in win_regions:
			for region in win_regions[id]:
				if region[0] > sp:
					new_seq += seq_db[id][sp: region[0]]
				new_seq += 'N'*(region[1]-region[0]+1)
				sp = region[1]+1
			if sp < len(seq_db[id]):
				new_seq += seq_db[id][sp: len(seq_db[id])]
		else:
			new_seq = seq_db[id]

		with open(id+'.tmp', 'w') as f_out:
			f_out.write(">%s\n%s\n"%(id, new_seq))


def quick_mask_genome(in_fa, in_bed, out_fa, th, ts):
	print("Reading fasta")
	seq_db, seq_id_list = read_fasta(in_fa)

	print("Reading bed")
	win_regions = merge_regions(read_bed(in_bed, th))

	task_per_thread = int(len(seq_id_list)/ts)
	
	task_list = []
	
	print("Masking genome")
	for i in range(0, ts):
		if i < ts-1:
			t = multiprocessing.Process(target=mask_fasta, args=(seq_id_list[i*task_per_thread: (i+1)*task_per_thread], seq_db, win_regions))
		else:
			t = multiprocessing.Process(target=mask_fasta, args=(seq_id_list[i*task_per_thread:], seq_db, win_regions))
		task_list.append(t)
	
	for t in task_list:
		t.start()
	
	for t in task_list:
		t.join()
	
	print("Merging")
	if os.path.exists(out_fa):
		os.remove(out_fa)
	for id in seq_id_list:
		os.system("cat "+id+".tmp >> "+out_fa)
		os.remove(id+".tmp")
	print("Success")


if __name__ == "__main__":
	if len(sys.argv) < 6:
		print("Usage: python "+sys.argv[0]+" <in_fasta> <in_bed> <out_fasta> <threshold> <threads>")
	else:
		in_fa = sys.argv[1]
		in_bed = sys.argv[2]
		out_fa = sys.argv[3]
		th = int(sys.argv[4])
		ts = int(sys.argv[5])
		quick_mask_genome(in_fa, in_bed, out_fa, th, ts)
