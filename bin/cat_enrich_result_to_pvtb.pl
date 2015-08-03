#!/usr/bin/perl -w 
#Author: Pengcheng Yang
#Email: yangpc@mail.biols.ac.cn

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname); 

my $usage= "\n$0 
--difcfg   \t<str>\tdiff enrich config
\t\tone file one line
--outDir   \t<str>\toutput directory

NOTE: this script is used to cat together all ther enrichment result into
one table file of the pvalue or the padjusted value
\n";

my($difcfg,$outDir);

GetOptions(
		"difcfg:s"=>\$difcfg,
		"outDir:s"=>\$outDir
);

die $usage if !defined $difcfg;

`mkdir -p $outDir` if ! -e $outDir;

my(@info,%pqv,%delid,%sample,%sets);
open I,$difcfg;
while(my $file=<I>){
	chomp $file;
	open I2,$file;
	#gSetID gSetDes Pvalue AdjustedPv x y n N EnrichDirect GeneIDs

	my $file_base=basename $file;
	$sample{$file_base}++;
	while(my $line1=<I2>){
		chomp $line1;
		@info=split /\t/,$line1;
		next if $info[0] =~ /gSetID/;
		next if $info[8] eq "Under";
		$pqv{$file_base}{$info[0]}{pv}=$info[2];
		$pqv{$file_base}{$info[0]}{qv}=$info[3];
		$delid{$info[0]}++ if $info[5] <5;
		$sets{$info[0]}=$info[1];
	}
	close I2;
}
close I;

my $difcfg_base = basename $difcfg;
open O,">","$outDir/$difcfg_base.pv.tb";
open O2,">","$outDir/$difcfg_base.qv.tb";
open O3,">","$outDir/$difcfg_base.id.title";
print O "Gsets\t",join("\t",sort keys %sample),"\n";
print O2 "Gsets\t",join("\t",sort keys %sample),"\n";
foreach my $setid (sort keys %sets){
	print O3 "$setid\t$sets{$setid}\n";
	next if exists $delid{$setid};
	print O "$setid";print O2 "$setid";
	foreach my $samp (sort keys %sample){
		if(exists $pqv{$samp}{$setid}){
			print O "\t$pqv{$samp}{$setid}{pv}";
			print O2 "\t$pqv{$samp}{$setid}{qv}";
		}else{
			print O "\t1";print O2 "\t1";
		}
	}
	print O "\n";print O2 "\n";
}
close O;close O2;


