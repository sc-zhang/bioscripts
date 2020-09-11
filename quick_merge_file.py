#!/usr/bin/env python
import sys, os
import multiprocessing


def sub_merge_file(f_list, s, e, tn):
	f_out = open('merge_tmp'+str(tn), 'w')
	for i in range(s, e):
		with open(f_list[i], 'r') as f_in:
			for line in f_in:
				f_out.write(line)
	f_out.close()


def merge_file(i_folder, f_t, m_f, tn):
	f_l = os.listdir(i_folder)
	merge_list = []
	for fn in f_l:
		if fn[-len(f_t):] == f_t:
			merge_list.append(os.path.join(i_folder, fn))
	task_per_prc = int(len(merge_list)/tn)
	task_list = []
	for i in range(0, tn):
		if i < tn-1:
			t = multiprocessing.Process(target=sub_merge_file, args=(merge_list, i*task_per_prc, (i+1)*task_per_prc, i))
		else:
			t = multiprocessing.Process(target=sub_merge_file, args=(merge_list, i*task_per_prc, len(merge_list), i))
		task_list.append(t)
	
	for t in task_list:
		t.start()

	for t in task_list:
		t.join()
	for i in range(0, tn):
		os.system("cat merge_tmp"+str(i)+" >> "+m_f)
		os.remove("merge_tmp"+str(i))


if __name__ == "__main__":
	if len(sys.argv) < 5:
		print("Usage: python "+sys.argv[0]+" <input_folder> <file_type> <merge_file> <threads>")
	else:
		prog, i_folder, f_t, m_f, t = sys.argv
		merge_file(i_folder, f_t, m_f, int(t))

