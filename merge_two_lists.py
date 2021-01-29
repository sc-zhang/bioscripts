#!/usr/bin/env python
import sys


def merge_two_lists(list1, col1, list2, col2, merged_list):
	list_db = {}
	with open(list2, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			list_db[data[col2]] = data
	
	with open(list1, 'r') as fin:
		with open(merged_list, 'w') as fout:
			for line in fin:
				data = line.strip().split()
				key = data[col1]
				if key in list_db:
					data.extend(list_db[key][:col2])
					data.extend(list_db[key][col2+1:])
				fout.write("%s\n"%('\t'.join(data)))


if __name__ == "__main__":
	if len(sys.argv) < 6:
		print("Usage: python %s <list1> <col1> <list2> <col2> <merged_list>"%sys.argv[0])
	else:
		list1, col1, list2, col2, merged_list = sys.argv[1:]
		col1 = int(col1)-1
		col2 = int(col2)-1
		merge_two_lists(list1, col1, list2, col2, merged_list)
