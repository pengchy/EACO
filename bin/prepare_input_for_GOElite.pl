#!/usr/bin/perl -w 
#Author: Pengcheng Yang
#Email: yangpc@mail.biols.ac.cn

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname); 

my $usage= "\n$0 
--idlist  \t<str>\tid list file, one file one line.
--outDir  \t<str>\toutput directory
--SysAbr  \t<str>\tSystem Code Abbreviation (default En)

\n";

my($idlist,$outDir,$SysAbr);

GetOptions(
		"idlist:s"=>\$idlist,
		"outDir:s"=>\$outDir,
		"SysAbr:s"=>\$SysAbr
);

die $usage if !defined $idlist;

$SysAbr ||= "En";

open I,$idlist;
while(my $file1=<I>){
	chomp $file1;
	my @info=split /\t/,$file1;
	my $file_base = basename $info[1];
	open O,">","$outDir/$file_base.txt";
	open I2,"$info[1]";
	while(<I2>){
		chomp;
		print O "$_\t$SysAbr\n";
	}
	close O;
	close I2;
}
close I;



