#!/usr/bin/env python
import sys


def rename_species(in_result, in_list, out_result):
	species_db = {}
	with open(in_list, 'r') as f_in:
		for line in f_in:
			data = line.strip().split()
			species_db[data[0]] = data[1]
	
	with open(in_result, 'r') as f_in:
		with open(out_result, 'w') as f_out:
			for line in f_in:
				data = line.strip().split()
				if data[0] == 'D':
					f_out.write(line)
				else:
					new_data = data[:8]
					f_out.write("\t".join(new_data))
					for i in range(8, len(data)):
						if data[i] in species_db:
							f_out.write("\t"+species_db[data[i]])
						else:
							f_out.write("\t"+data[i])
					f_out.write("\n")


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Notice: rename species in result.ErrorCorr.txt with list file")
		print("Usage: python "+sys.argv[0]+" <in_result> <in_list> <out_result>")
	else:
		in_result = sys.argv[1]
		in_list = sys.argv[2]
		out_result = sys.argv[3]
		rename_species(in_result, in_list, out_result)

