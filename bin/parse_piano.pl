#!/usr/bin/perl -w 
#Author: Pengcheng Yang
#Email: yangpc@mail.biols.ac.cn

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname); 
use GACP qw(parse_config);

my $usage= "\n$0 
--bincfg  \t<str>\tbin configure file, in the format: Name = file_path
--outDir  \t<str>\toutput directory
--pianof  \t<str>\tpiano file in the format: sampleID<tb>piano_file

\n";

my($bincfg,$formt,$outDir,$pianof);

GetOptions(
		"bincfg:s"=>\$bincfg,
		"outDir:s"=>\$outDir,
		"pianof:s"=>\$pianof
);

die $usage if !defined $pianof;

my $pianof_base=basename $pianof;
my(@info,$file,%pv,%padj,%pos,$pos_up,$pos_dn,$pos_up_adj,$pos_dn_adj,%gsets);
open I,$pianof;
while(my $lin1=<I>){
	chomp $lin1;
	@info=split /\t/,$lin1;
	$info[0] =~ s/\-/_/g;
	$file = basename $info[1];
	open I2,$info[1];
	while(my $lin2=<I2>){
		next if $lin2 !~ /\w/;
		my @info2=split /\t/, $lin2;
		if($info2[0] eq "Name"){
			my $i;
			%pos=map{$_=>$i++}@info2;
			$pos_up_adj=$pos{"p adj (dist.dir.up)"};
			$pos_dn_adj=$pos{"p adj (dist.dir.dn)"};
			$pos_up=$pos{"p (dist.dir.up)"};
			$pos_dn=$pos{"p (dist.dir.dn)"};
			next;
		}
		$pv{$info[0]}{$info2[0]}=$info2[$pos_up]."\t".$info2[$pos_dn];
		$padj{$info[0]}{$info2[0]}=$info2[$pos_up_adj]."\t".$info2[$pos_dn_adj];
		$gsets{$info2[0]}++;
#print "$lin1\n"; print "$lin2\n"; print $pv{$info[0]}{$info2[0]},"\n"; last;
	}
	close I2;
#last;
}
close I;

open O,">","$outDir/$pianof_base.pv.tb";
open O2,">","$outDir/$pianof_base.padj.tb";
my @files = sort keys %pv;
print O "Gsets";
print O2 "Gsets";
map{print O "\t$_.Up\t$_.Down" }@files;
map{print O2 "\t$_.Up\t$_.Down" }@files;
print O "\n";
print O2 "\n";

for my $gset (sort keys %gsets){
	print O "$gset";
	print O2 "$gset";
	for $file (@files){
		if(exists $pv{$file}{$gset}){
			print O "\t$pv{$file}{$gset}";
			print O2 "\t$padj{$file}{$gset}";
		}else{
			print O "\tNA\tNA";
			print O2 "\tNA\tNA";
		}
	}
	print O "\n";
	print O2 "\n";
}
close O;close O2;

