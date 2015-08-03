#!/usr/bin/perl -w 
#Author: Pengcheng Yang
#Email: yangpc@mail.biols.ac.cn

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname); 

my $usage= "\n$0 
--gmt     \t<str>\tgmt file
--outF    \t<str>\ttwo column file

NOTE: convert file in gmt format into two column format

\n";

my($gmt,$outF);

GetOptions(
		"gmt:s"=>\$gmt,
		"outF:s"=>\$outF
);

die $usage if !defined $gmt;

my(@info,%data);
open I,$gmt;
open O,">","$outF";
while(<I>){
	chomp;
	@info=split /\t/;
	map{print O "$_\t$info[0]\n";}@info[2..$#info];
}
close O;
close I;



