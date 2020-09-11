#!/usr/bin/env python
import sys
import multiprocessing


def write_cmd(fn, cmd_list):
	print("\tWriting %s"%fn)
	with open(fn, 'w') as fout:
		fout.write("".join(cmd_list))


def split_cmd(in_cmd, np, out_str, ts):
	print("Loading cmds")
	cmd_list = []
	with open(in_cmd, 'r') as fin:
		for line in fin:
			cmd_list.append(line)
	
	print("Splitting commands")
	pool = multiprocessing.Pool(processes=ts)
	cmd_per_file = int(round(len(cmd_list)/np, 0))
	for i in range(0, np):
		fn = out_str%(i+1)
		if i < np-1:
			pool.apply_async(write_cmd, (fn, cmd_list[i*cmd_per_file: (i+1)*cmd_per_file],))
		else:
			pool.apply_async(write_cmd, (fn, cmd_list[i*cmd_per_file:],))
	pool.close()
	pool.join()
	print("Finished")


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python "+sys.argv[0]+" <in_cmd_file> <num_parts> <out_str> <threads>")
		print("\t<out_str> is a string contain %d as file index, like run_%d.sh")
	else:
		in_cmd, np, out_str, ts = sys.argv[1:]
		split_cmd(in_cmd, int(np), out_str, int(ts))

