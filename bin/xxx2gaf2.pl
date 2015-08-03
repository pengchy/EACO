#!/usr/bin/perl -w 
#Author: Pengcheng Yang
#Email: yangpc@mail.biols.ac.cn

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname); 
use GACP qw(parse_config);

my $usage= "\n$0 
--bincfg  \t<str>\tbin configure file, in the format: Name = file_path
--xxx     \t<str>\tinput file
--inf     \t<str>\tinput file format
\t\ttb2  two column, geneID<tb>GOid (default)
\t\tmtb  multiple column, geneID<tb>GOid<tb>GOid (unimplemented)
\t\tgmt  gmt format, GOid<tb>GOdescription<tb>geneID<tb>geneID<tb>...
--anno    \t<str>\tgene annotation information. geneID<tb>geneDescription
--Onto    \t<str>\tgene ontology table file. GOid<tb>Class<tb>Desc
--outf    \t<str>\toutput file
--tax     \t<str>\tNCBI taxonomy ID (default 7004, which is locusta migratoria)

NOTE: this script is used to produce gaf2 file format from various file
the gaf2 file format specification can be get from:
http://www.geneontology.org/page/go-annotation-file-gaf-format-20

Column	Content	Required?	Cardinality	Example
1	DB	required	1	UniProtKB
2	DB Object ID	required	1	P12345
3	DB Object Symbol	required	1	PHO3
4	Qualifier	optional	0 or greater	NOT
5	GO ID	required	1	GO:0003993
6	DB:Reference (|DB:Reference)	required	1 or greater	PMID:2676709
7	Evidence Code	required	1	IMP
8	With (or) From	optional	0 or greater	GO:0000346
9	Aspect	required	1	F
10	DB Object Name	optional	0 or 1	Toll-like receptor 4
11	DB Object Synonym (|Synonym)	optional	0 or greater	hToll|Tollbooth
12	DB Object Type	required	1	protein
13	Taxon(|taxon)	required	1 or 2	taxon:9606
14	Date	required	1	20090118
15	Assigned By	required	1	SGD
16	Annotation Extension	optional	0 or greater	part_of(CL:0000576)
17	Gene Product Form ID	optional	0 or 1	UniProtKB:P12345-2

\n";

#get the parameters and predefine variables
my($bincfg,$xxx,$inf,$db,$outf,$anno,$Onto,$tax);
GetOptions(
		"bincfg:s"=>\$bincfg,
		"xxx:s"=>\$xxx,
		"inf:s"=>\$inf,
		"outf:s"=>\$outf,
		"anno:s"=>\$anno,
		"Onto:s"=>\$Onto,
		"tax:s"=>\$tax
);
die $usage if !defined $xxx;
$inf ||= "tb2";
$tax ||= "7004";

my(@info,%gen2go);
my($gid,$goid,%gaf,%an,%On);

#read in annotation and Ontology
open I,$anno;
while(<I>){
	chomp;
	@info=split /\t/;
	$an{$info[0]}{desc}=$info[1];
	$an{$info[0]}{syn}=$info[2];
	if($info[1] eq "NA"){
		$an{$info[0]}{desc}="hypothetical protein";
	}
	if($info[2] eq "NA"){
		$an{$info[0]}{syn}=$info[0];
	}
}
close I;

open I,$Onto;
while(<I>){
	chomp;
	@info=split /\t/;
	$On{$info[0]}{class}="P" if $info[1] eq "BP";
	$On{$info[0]}{class}="F" if $info[1] eq "MF";
	$On{$info[0]}{class}="C" if $info[1] eq "CC";
	$On{$info[0]}{def}=$info[2];
}
close I;

##predefine gaf
$gaf{1}="LocustGenomeProject";
$gaf{4}=0;
$gaf{6}="NA";
$gaf{7}="IEA";
$gaf{8}="NA";
$gaf{9}="NA";
$gaf{11}="NA";
$gaf{12}="gene";
$gaf{13}= "taxon:$tax";
my $dat=`date +"%Y%m%d"`;chomp $dat;
$gaf{14}=$dat;
$gaf{15}="LocustGenomeProject";
if($inf eq "tb2"){
	open I,$xxx;
	while(<I>){
		chomp;
		@info=split /\t/;
		$gen2go{$info[0]}{$info[1]}++;
	}
	close I;
	open O,">","$outf";
	foreach $gid (sort keys %gen2go){
		$gaf{2}=$gid;$gaf{3}=$an{$gid}{syn};
		$gaf{10}=$an{$gid}{desc};
		foreach $goid (sort keys %{$gen2go{$gid}}){
			$gaf{5}=$goid;
			$gaf{9}=$On{$goid}{class} if exists $On{$goid}{class};
			print O $gaf{1};
			map{print O "\t",$gaf{$_} }(2..15);
			print O "\n";
		}
	}
}
close O;
