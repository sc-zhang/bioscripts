#!/usr/bin/python
import os, sys
import ConfigParser


def get_config(f_conf, section, key):
	config = ConfigParser.ConfigParser()
	config.read(f_conf)
	return config.get(section, key)


def get_options(Argvs):
	i = 0
	opts = {}
	fasta_file = ''
	gff_files = []
	bed_files = []
	bam_files = []
	bw_files = []
	vcf_files = []
	conf_file = ''
	chrs_file = ''
	tn = '1'
	while i<len(Argvs):
		if Argvs[i] == '-f':
			fasta_file = Argvs[i+1]
		elif Argvs[i] == '--gff':
			gff_files.append(Argvs[i+1])
		elif Argvs[i] == '--bed':
			bed_files.append(Argvs[i+1])
		elif Argvs[i] == '--bam':
			bam_files.append(Argvs[i+1])
		elif Argvs[i] == '--bw':
			bw_files.append(Argvs[i+1])
		elif Argvs[i] == '--conf':
			conf_file = Argvs[i+1]
		elif Argvs[i] == '--bam2bw':
			chrs_file = Argvs[i+1]
		elif Argvs[i] == '-t':
			tn = Argvs[i+1]
		elif Argvs[i] == '--vcf':
			vcf_files.append(Argvs[i+1])
		elif Argvs[i] == '-h' or Argvs[i] == '--help':
			print("Usage: python "+sys.argv[0]+" -f <fasta_file> [--gff <gff_file> --bed <bed_file> --bam <bam_file> --bw <bigwig_file> --conf <config_file>]")
			exit(0)
		i += 2
	opts['fasta_file'] = fasta_file
	opts['gff_files'] = gff_files
	opts['bed_files'] = bed_files
	opts['bam_files'] = bam_files
	opts['vcf_files'] = vcf_files
	opts['bw_files'] = bw_files
	opts['bam2bw'] = chrs_file
	opts['conf_file'] = conf_file
	opts['tn'] = tn
	return opts


def simple_jbrowser(opts):
	samtools = ''
	jbrowser = ''
	bam2wig = ''
	if opts['conf_file'] != '':
		conf_file = opts['conf_file']
		if os.path.isfile(conf_file):
			samtools = get_config(conf_file, "path", "samtools")
			bam2wig = get_config(conf_file, "path", "bam2wig")
			jbrowser = get_config(conf_file, "path", "JBrowser")
			wig2bw = get_config(conf_file, "path", "wig2bw")
	else:
		print("No configure file")
		exit(0)
	
	if samtools != '' and samtools[-1] != '/':
		samtools += '/'
	if jbrowser != '' and jbrowser[-1] != '/':
		jbrowser += '/'
	if bam2wig != '' and bam2wig[-1] != '/':
		bam2wig += '/'
	if wig2bw != '':
		os.system("export PATH="+wig2bw+":$PATH")

	print("Preparing reference sequences")
	if opts['fasta_file'] == '':
		print("No reference sequences")
		exit(0)
	os.system(jbrowser+"prepare-refseqs.pl --fasta "+opts['fasta_file'])

	print("Preparing gffs")
	for gff in opts['gff_files']:
		if gff != '':
			os.system(jbrowser+"flatfile-to-json.pl --gff "+gff+" --trackType CanvasFeatures --trackLabel "+gff.split('.')[0])
	
	print("Preparing beds")
	for bed in opts['bed_files']:
		if bed != '':
			os.system(jbrowser+"flatfile-to-json.pl --bed "+bed+" --trackType CanvasFeatures --trackLabel "+bed)
	
	print("Preparing vcf")
	for vcf in opts['vcf_files']:
		if vcf[-3:].lower() != '.gz':
			os.system("bgzip "+vcf)
			vcf = vcf+".gz"
		if os.path.exists(vcf+".tbi") == False:
			os.system("tabix -p vcf "+vcf)
		with open("data/tracks.conf", "a") as f_track:
			f_track.write("[tracks."+vcf.replace('.', '_')+"]\nstoreClass = JBrowse/Store/SeqFeature/VCFTabix\nurlTemplate = ../"+vcf+"\ncategory = VCF\ntype = JBrowse/View/Track/CanvasVariants\nkey  = "+vcf.replace('.', '_')+"\n")

	print("Preparing bam")
	for bam in opts['bam_files']:
		if not bam.endswith('sorted.bam'):
			sorted_bam = bam+".sorted.bam"
			indexed_bam = sorted_bam+".bai"
			if not os.path.exists(sorted_bam):
				os.system(samtools+"samtools sort -@ "+opts['tn']+" -o "+sorted_bam+" "+bam)
			if not os.path.exists(indexed_bam):
				os.system(samtools+"samtools index "+sorted_bam)
		else:
			sorted_bam = bam
			indexed_bam  = sorted_bam+".bai"
			if not os.path.exists(indexed_bam):
				os.system(samtools+"samtools index "+sorted_bam)
		with open("data/tracks.conf", "a") as f_track:
			f_track.write("[tracks."+bam.replace('.', '_')+"]\nstoreClass = JBrowse/Store/SeqFeature/BAM\nurlTemplate = ../"+sorted_bam+"\nbaiUrlTemplate = ../"+indexed_bam+"\ncategory = NGS\ntype = JBrowse/View/Track/Alignments2\nkey = "+bam.replace('.', '_')+"\n")
		if opts['bam2bw'] != '':
			os.system(bam2wig+"bam2wig.py -i "+sorted_bam+" -s "+opts['bam2bw']+" -o "+bam)
			opts['bw_files'].append(bam+".bw")
			os.remove(bam+".wig")
	
	print("Preparing bigwig")
	for bw in opts['bw_files']:
		with open("data/tracks.conf", "a") as f_track:
			f_track.write("[tracks."+bw.replace('.', '_')+"]\nstoreClass = JBrowse/Store/SeqFeature/BigWig\nurlTemplate = ../"+bw+"\ncategory = Quantitative\ntype = JBrowse/View/Track/Wiggle/XYPlot\nkey = "+bw.replace('.', '_')+"\n")
	print("Finished")


if __name__ == "__main__":
	if len(sys.argv) == 1:
		print("Usage: python "+sys.argv[0]+" -f <ref_fasta> --conf <config_file> [--gff <gff_file> --bed <bed_file> --vcf <vcf_file> --bam <bam_file> --bw <bigwig_file> --bam2bw <chrs_length_file> -t <threads>]")
	else:
		opts = get_options(sys.argv[1:])
		simple_jbrowser(opts)
