#!/usr/bin/perl -w 
#Author: Pengcheng Yang
#Email: yangpc@mail.biols.ac.cn

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname); 

my $usage= "\n$0 
--dt      \t<str>\tthree column data file  Set1<tb>Set2<tb>Value
--mt      \t<str>\toutput matrix

NOTE: used to convert the three column data into matrix format

\n";

my($dt,$mt);

GetOptions(
		"dt:s"=>\$dt,
		"mt:s"=>\$mt
);

die $usage if !defined $dt;

my(@info,%mtr);
open I,"$dt";
while(<I>){
	chomp;
	@info=split /\t/;
	if($info[2] < 0.0000001) {$info[2] = 0;}
	$mtr{$info[0]}{$info[1]}=$info[2];
	$mtr{$info[1]}{$info[0]}=$info[2];
}
close I;

open O,">","$mt";
my @ids = sort {$a cmp $b} (keys %mtr);
print O join("\t",@ids),"\n";
foreach my $id1 (@ids){
	print O "$id1";
	foreach my $id2 (@ids){
		print O "\t$mtr{$id1}{$id2}";
	}
	print O "\n";
}
close O;

