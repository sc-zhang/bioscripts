#!/usr/bin/perl -w

use Getopt::Std;
getopts "i:q:e:";


if ((!defined $opt_i)|| (!defined $opt_e) || (!defined $opt_q)) {
    die "************************************************************************
    Usage: perl $0 -i contig.fasta -q fastq_dir -e enzyme 
      -h : help and usage.
      -q : fastq_dir 
      -i : contig.fasta
      -e : HINDIII/MBOI
      -s : step, int (from 1 to X)
************************************************************************\n";
}

system("rm -rf scripts");
system("mkdir scripts");

my $refSeq = $opt_i;
my $enzy   = $opt_e;
system("ln -s $refSeq ./draft.asm.fasta");
system("bwa index -a bwtsw draft.asm.fasta");
system("samtools faidx draft.asm.fasta");

while(my $file = glob "$opt_q/*R1_*\.*"){
	my ($pre, $id, $suf) = ($1, $2, $3) if($file=~/\/(\S+)_R1_(\d+)\.(\S+)/);
	my $r1     = $opt_q."/".$pre."_R1_".$id.".".$suf;
	my $r2     = $opt_q."/".$pre."_R2_".$id.".".$suf;
	my $outfile = "run_".$id.".sh";
	my $sai1   = "sample".$id."_R1.sai";
	my $sai2   = "sample".$id."_R2.sai";
	my $bwa_aln = "sample".$id.".bwa_aln.sam";
	my $PE_bam = "sample".$id.".bwa_aln.REduced.paired_only.bam";
	my $filter_sam = "filter".$id.".sam";
	open(my $out, "> scripts/$outfile") or die"";
	print $out "#!/bin/bash\n";
	print $out "bwa aln -t 6 draft.asm.fasta $r1 >$sai1\n";
	print $out "bwa aln -t 6 draft.asm.fasta $r2 > $sai2\n";
	print $out "bwa sampe draft.asm.fasta $sai1 $sai2 $r1 $r2 > $bwa_aln\n";
	print $out "PreprocessSAMs.pl $bwa_aln draft.asm.fasta $enzy\n";
	print $out "filterBAM_forHiC.pl $PE_bam $filter_sam\n";
	close $out;
		system("qsub -S /bin/bash -cwd -pe mpi 10 -q all.q scripts/$outfile");
	}







