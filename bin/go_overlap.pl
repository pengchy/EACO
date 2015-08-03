#!/usr/bin/perl -w 
#Author: Pengcheng Yang
#Email: yangpc@mail.biols.ac.cn

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname); 
use GACP qw(parse_config);

my $usage= "\n$0 
--bincfg  \t<str>\tbin configure file, in the format: Name = file_path
--go      \t<str>\tgene association file in the format geneID<tb>GOid
--outf    \t<str>\toutput file

\n";

my($bincfg,$go,$outf);

GetOptions(
		"bincfg:s"=>\$bincfg,
		"go:s"=>\$go,
		"outf:s"=>\$outf
);

die $usage if !defined $go;

my(@info,%go2gen);
open I,$go;
while(<I>){
	chomp;
	@info=split /\t/;
	$go2gen{$info[1]}{$info[0]}++;
}
close I;

open O,">","$outf";
print O "GO1\tGO2\tGO1Num\tGO2Num\tOverlap\t";
my @gos = sort keys %go2gen;
for(my $i=0;$i<$#gos;$i++){
	for(my $j=$i+1;$j<=$#gos;$j++){
		print O ""
	}
}

