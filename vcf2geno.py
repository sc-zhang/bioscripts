#!/usr/bin/env python
import sys


def help_message():
	print("Notice: convert vcf file to geno file for ABBABABAwindows.py")
	print("Usage: python "+sys.argv[0]+" -i <input_vcf> -o <output_vcf> -q/--quality <min_qual> -f/--filter <filter_type> <min_value>")


def Getopts(ARGV):
	opts = ARGV
	input_file = ''
	output_file = ''
	quality = 0
	filter_db = {}
	for i in range(0, len(opts)):
		if opts[i] == '-i':
			input_file = opts[i+1]
		if opts[i] == '-o':
			output_file = opts[i+1]
		if opts[i] == '-q' or opts[i] == '--quality':
			quality = int(opts[i+1])
		if opts[i] == '-f' or opts[i] == '--filters':
			j = i+1
			while(j<len(opts) and opts[j][0] != '-'):
				filter_db[opts[j]] = int(opts[j+1])
				j += 2
		if opts[i] == '-h' or opts[i] == '--help':
			help_message()
			exit(0)
	if input_file == '' or output_file == '':
		help_message()
		exit(0)
	return input_file, output_file, quality, filter_db


def vcf_parser(input_file, output_file, quality, filter_db):
	header = []
	with open(input_file, 'r') as f_in:
		with open(output_file, 'w') as f_out:
			for line in f_in:
				data = line.strip().split()
				if line[:2]  == '##':
					continue
				if line[0] == '#':
					header = data
					f_out.write(header[0]+"\t"+header[1])
					for i in range(9, len(header)):
						f_out.write("\t"+header[i])
					f_out.write('\n')
				else:
					chrn = data[0]
					pos = data[1]
					ref = data[3]
					alt = data[4]
					qual = float(data[5])
					if qual < quality:
						continue
					info = data[7].split(';')
					is_skip = False
					for item in info:
						type = item.split('=')[0]
						value = item.split('=')[1]
						if value.isdigit():
							value = float(value)
						else:
							continue
						if type in filter_db:
							if value < filter_db[type]:
								is_skip = True
								break
					if is_skip:
						continue
					f_out.write(chrn+"\t"+pos)
					for i in range(9, len(data)):
						species = data[i].split(':')[0]
						species = species.replace('0', ref)
						species = species.replace('1', alt)
						species = species.replace('.', 'N')
						f_out.write("\t"+species)
					f_out.write('\n')



def main():
	input_file, output_file, quality, filter_db = Getopts(sys.argv[1:])
	vcf_parser(input_file, output_file, quality, filter_db)


if __name__ == "__main__":
	if len(sys.argv) == 1:
		help_message()
	else:
		main()
