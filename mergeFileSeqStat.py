#!/usr/bin/env python
import sys
import os
import pysam
import gzip


def read_fasta(in_file):
	if in_file.split('.')[-1] == 'gz':
		fin = gzip.open(in_file, 'rt')
	else:
		fin = open(in_file, 'r')
	seq_len = []
	for line in fin:
		if line[0] == '>':
			seq_len.append(0)
		else:
			seq_len[-1] += len(line.strip())
	fin.close()
	return seq_len


def read_fastq(in_file):
	if in_file.split('.')[-1] == 'gz':
		fin = gzip.open(in_file, 'rt')
	else:
		fin = open(in_file)
	seq_len = []
	cnt = 0
	for line in fin:
		if cnt%4 == 1:
			seq_len.append(len(line.strip()))
		cnt += 1
	fin.close()
	return seq_len


def read_bam(in_file):
	seq_len = []
	'''
	with os.popen("samtools view %s"%in_file, 'r') as fin:
		for line in fin:
			seq_len.append(len(line.strip().split()[9]))
	'''
	with pysam.AlignmentFile(in_file, 'rb', check_sq=False) as fin:
		for line in fin:
			seq_len.append(line.query_length)
	return seq_len


def check_file_type(in_file):
	data = in_file.split('.')
	if data[-1] == 'gz':
		with gzip.open(in_file, 'rt') as fin:
			for line in fin:
				break
			if line[0] == '>':
				return "fa"
			elif line[0] == '@':
				return "fq"
			else:
				return ""
	elif data[-1] == 'bam':
		return "bam"
	else:
		with open(in_file, 'r') as fin:
			for line in fin:
				break
			if line[0] == '>':
				return "fa"
			elif line[0] == '@':
				return "fq"
			else:
				return ""


def seq_stat(in_files, out_stat):
	seq_len = []
	for in_file in in_files.split(','):
		print("Reading %s"%in_file)
		file_type = check_file_type(in_file)
		if file_type == 'fa':
			seq_len.extend(read_fasta(in_file))
		elif file_type == 'fq':
			seq_len.extend(read_fastq(in_file))
		elif file_type == 'bam':
			seq_len.extend(read_bam(in_file))
		else:
			print("Unsupport file type")
			continue
	seq_len = sorted(seq_len, reverse=True)
	seq_cnt = len(seq_len)
	if seq_cnt == 0:
		print("No seq load")
		sys.exit()
	min_len = seq_len[-1]
	max_len = seq_len[0]
	total_size = sum(seq_len)
	ave_len = total_size*1.0/seq_cnt
	n_threshold = []
	n_values = []
	n_labels = []
	for i in range(90, 40, -10):
		n_threshold.append(i/100.0*total_size)
		n_values.append(0)
		n_labels.append("N%d:\t"%i)
	cur_size = 0
	cnt_500 = 0
	cnt_2k = 0
	for i in range(0, seq_cnt):
		cur_size += seq_len[i]
		for j in range(0, len(n_values)):
			if n_values[j] == 0 and cur_size >= n_threshold[j]:
				n_values[j] = seq_len[i]
		if seq_len[i] > 500:
			cnt_500 += 1
		if seq_len[i] > 2000:
			cnt_2k += 1
	n_info = []
	for i in range(0, len(n_values)):
		n_info.append("%s%d\n"%(n_labels[i], n_values[i]))
	info = "number of seq:\t%d\nmin length:\t%d\nmax length:\t%d\ntotal size:\t%d\n%sAverage length:\t%d\nTotal number (>500bp):\t%d\nTotal number (>2000bp):\t%d"%(seq_cnt, min_len, max_len, total_size, ''.join(n_info), ave_len, cnt_500, cnt_2k)
	if out_stat == "":
		print(info)
	else:
		with open(out_stat, 'w') as fout:
			fout.write("%s\n"%info)


if __name__ == "__main__":
	if len(sys.argv) < 2:
		print("Usage: python %s <in_files> [out_stat]"%sys.argv[0])
	else:
		if len(sys.argv) == 2:
			in_files = sys.argv[1]
			out_stat = ""
		else:
			in_files, out_stat = sys.argv[1:]
		seq_stat(in_files, out_stat)

