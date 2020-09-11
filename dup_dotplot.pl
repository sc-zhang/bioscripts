#!/usr/bin/perl -w
use Getopt::Std;
getopts "g:r:q:n:t:";

if ((!defined $opt_g)or(!defined $opt_r)or(!defined $opt_q)or(!defined $opt_n)){
	die"******************************************************************************************
    Usage: perl $0 -g reference_genome -r ref_id -q query_id -n number_of_dup -t threads
        ref_id: reference cds and bed name, like: Sb, Sb.cds and Sb.bed must exist
        query_id: query cds and bed name, like: Os
        number_of_dup: number of duplications
        threads: default 1
******************************************************************************************\n";
}

my $genome = $opt_g;
my $ref_name = $opt_r;
my $qry_name = $opt_q;
my $dup_n = $opt_n;
if (!defined $opt_t){
	$threads = "1";
}
else{
	$threads = $opt_t;
}

my %sbcdsdb;
open(IN, $ref_name.".cds") or die"";
while(<IN>){
	chomp;
	if(/>/){
		$gene = $_;
		$gene =~ s/\s.*//g;
		$gene =~ s/>//g;
	}else{
		$sbcdsdb{$gene} .= $_;
		}
	}
close IN;

system("gmap_build -D . -d DB ".$genome);
system("gmap -D . -d DB -t ".$threads." -f 2 -n ".$dup_n." ".$ref_name.".cds"." > ".$qry_name.".gff3");

open(CDS, "> ".$qry_name.".cds") or die"";
open(BED, "> ".$qry_name.".bed") or die"";
my $count = 0;
open(IN, "grep 'gene' ".$qry_name.".gff3|") or die"";
while(<IN>){
	chomp;
	$count++;
	my @data = split(/\s+/,$_);
	my $gid  = $data[0];
	my $a    = $data[3];
	my $b    = $data[4];
	my $gene = $1 if(/Name=(\S+)/);
	   $gene =~ s/;.*//g;
	my $cds  = $sbcdsdb{$gene};
	   $gene = $gene."_".$count;
	print CDS ">$gene\n$cds\n";
	print BED "$gid	$a	$b	$gene	0	$data[6]\n";
	}
close IN;
close CDS;
close BED;

system(" python -m jcvi.compara.catalog ortholog ".$qry_name." ".$ref_name);

