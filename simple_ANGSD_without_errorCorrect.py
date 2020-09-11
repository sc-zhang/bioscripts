#!/usr/bin/env python
import sys, os


def help_message():
	print("Usage: python "+sys.argv[0]+" -l <species.list> -r <region> [-out <out_group_name> -p <bam_path>]")


def parse_options(ARGV):
	opt_dict = {}
	for i in range(0, len(ARGV), 2):
		opt = ARGV[i]
		if opt == '-l':
			opt_dict['list'] = ARGV[i+1]
		elif opt == '-p':
			opt_dict['path'] = ARGV[i+1]
		elif opt == '-r':
			opt_dict['regions'] = ARGV[i+1]
		elif opt == '-out':
			opt_dict['out_name'] = ARGV[i+1]
		elif opt == '-h':
			help_message()
			exit(0)
	return opt_dict


def run_abbababa(opts):
	if 'path' not in opts:
		bam_path = './'
	else:
		bam_path = opts['path']
	bam_files = []
	for fn in os.listdir(bam_path):
		if fn[-4:] == '.bam':
			bam_files.append(os.path.join(bam_path, fn))
	
	print("Indexing bams")
	for fn in bam_files:
		fn = fn.split('/')[-1]
		if os.path.isfile(os.path.join('./', fn+'.bai')) == False:
			cmd = 'samtools index '+fn
			print("Running command: "+cmd)
			os.system(cmd)
	
	print("Done\nReading list")
	list_db = {}
	with open(opts['list'], 'r') as f_in:
		for line in f_in:
			data = line.strip().split()
			if data[1] not in list_db:
				list_db[data[1]] = {}
				list_db[data[1]]['name'] = []
				list_db[data[1]]['path'] = []
			for fn in bam_files:
				if data[0] in fn:
					list_db[data[1]]['path'].append(fn)
					list_db[data[1]]['name'].append(data[0])
	
	if 'out_name' not in opts:
		out_name = "Outgroup"
	else:
		out_name = opts['out_name']
	
	print("Done\nGenerate bam.filelist sizeFile.size popNames.name bamWithErrors.filelist errorList.error")
	with open("bam.filelist", "w") as f_list:
		with open("sizeFile.size", "w") as f_size:
			with open("popNames.name", "w") as f_pop:
				for subgroup in list_db:
					if subgroup != out_name:
						f_pop.write(subgroup+'\n')
						f_list.write('\n'.join(list_db[subgroup]['path'])+'\n')
						group_size = len(list_db[subgroup]['name'])
						f_size.write(str(group_size)+'\n')
				f_list.write(list_db[out_name]['path'][0]+'\n')
				f_pop.write(out_name+'\n')
				f_size.write('1\n')
	
	print("Done\nDo abbababa")
	if "regions" not in opts:
		regions_file = "regions.txt"
	else:
		regions_file = opts["regions"]

	cmd = "ANGSD -doAbbababa2 1 -bam bam.filelist -sizeFile sizeFile.size -doCounts 1 -out bam.Angsd -rf "+regions_file+" -useLast 1 -minQ 20 -minMapQ 30"
	print("Running command: "+cmd)
	os.system(cmd)
	
	print("Done\n")
	
	cmd = "Rscript /public1/home/zsc/software/angsd/DSTAT angsdFile=\"bam.Angsd\" out=\"result\" sizeFile=sizeFile.size nameFile=popNames.name"
	print("Running command: "+cmd)
	os.system(cmd)
	print("Done\nSuccess")


if __name__ == "__main__":
	if len(sys.argv) == 1:
		help_message()
	else:
		opts = parse_options(sys.argv[1:])
		run_abbababa(opts)
