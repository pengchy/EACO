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
--glist   \t<str>\tgene list file, one line one file
--gsets   \t<str>\tgene set file, in the format: SetName<tb>SetFile
--method  \t<str>\talgorithm used to calculate the enrichment score, currently support:
\t\tfisher, stouffer, reporter, tailStrength, wilcoxon, mean, median, sum, maxmean, gsea or page.
\t\tevery method can be concatenated by comma. default is gsea

NOTE: this script is used to run piano in batch mode
\n";

my($bincfg,$outDir,$glist,$gsets,$method);

GetOptions(
		"bincfg:s"=>\$bincfg,
		"glist:s"=>\$glist,
		"gsets:s"=>\$gsets,
		"outDir:s"=>\$outDir,
		"method:s"=>\$method
);

die $usage if !defined $glist;

my $Bin=dirname $0;
my $run_piano = "perl $Bin/run_piano.pl";

$method ||= "gsea";

my(@info,%sets);
open I,$gsets;
while(<I>){
	chomp;
	@info=split /\t/;
	$sets{$info[0]}=$info[1];
}
close I;

my($id);
`mkdir -p $outDir` if ! -e $outDir;
open I1,$glist;
open O,">","$outDir/run_piano.sh";
open O2,">","$outDir/parse_piano.sh";
print O2 "#cat together the piano results\n";
while(<I1>){
	chomp;
	@info=split /\t/;
	foreach my $set (sort keys %sets){
		print O "$run_piano --outDir $info[1].enrich --glevel $info[1] --gset $sets{$set} --SetName $set --method $method\n";
	}
	foreach $id (split /,/,$method){
		print O2 "cat $info[1].enrich/$info[0].piano.*.$id |head -n 1 > $outDir/$info[0].$id\n";
		print O2 qq#cat $info[1].enrich/$info[0].piano.*.$id |grep -v "^Name"  >> $outDir/$info[0].$id\n#;
	}
}
close I1;
close O;
close O2;

open O2,">>","$outDir/parse_piano.sh";
foreach $id (split /,/,$method){
	open I1,$glist;
	open O,">","$outDir/enrich.$id";
	while(<I1>){
		chomp;@info=split /\t/;
		print O "$info[0]\t$outDir/$info[0].$id\n";
	}
	close I1;
	close O;

	print O2 "\n#parse piano file\n";
	print O2 "perl $Bin/parse_piano.pl --outDir $outDir/ --pianof $outDir/enrich.$id\n";
}
close O2;


