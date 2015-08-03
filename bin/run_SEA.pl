#!/usr/bin/perl -w 
#Author: Pengcheng Yang
#Email: yangpc@mail.biols.ac.cn

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname); 

my $usage= qq#\n$0 
--gmt     \t<str>\tgene sets file in gmt format
--outDir  \t<str>\toutput directory
--supplyID\t<str>\tids to enrich, one id one line
--method  \t<str>\tcurrently support
\t\tFisherChiSquare, HyperGeometric (default FisherChiSquare)
--bakid   \t<str>\tbackground id file, one id one line, if undefined the all
\t\tids in the gmt will be used as background.

--padjM   \t<str>p adjust method: "holm","hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"
\t\t(default fdr)

NOTE: this script is used to run enrichment analysis for the supplied id list
using the SEA method.
The SEA stands for Singular Enrichment Analysis.

\n#;

my($gmt,$outDir,$supplyID,$method,$bakid,$padjM);

GetOptions(
		"gmt:s"=>\$gmt,
		"outDir:s"=>\$outDir,
		"supplyID:s"=>\$supplyID,
		"method:s"=>\$method,
		"bakid:s"=>\$bakid,
		"padjM:s"=>\$padjM
);

die $usage if !defined $gmt;

my $id_base = basename $supplyID;
my $dir = dirname $0;
my $rscp;

`mkdir -p $outDir` if ! -e "$outDir";

$bakid ||= "NullFile";
$padjM ||= "fdr";
$method ||= "FisherChiSquare";


$rscp=qq#
date()
source("$dir/EnrichSGA.R")
source("$dir/sort.data.frame.R")
source("$dir/Fisher.Chi.test.R")
supplyID <- scan("$supplyID",what="character",sep="\\n")
if(identical("$bakid","NullFile")){
	res.dt <- EnrichSGA(gmt="$gmt",supplyID=supplyID,p.adjust.methods="$padjM",
			test.method="FisherChiSquare",enrichFile="$outDir/$id_base.difsea")
}else{
	univerID <- scan("$bakid",what="character",sep="\\n")
	res.dt <- EnrichSGA(gmt="$gmt",supplyID=supplyID,univerID=univerID,p.adjust.methods="$padjM",
			test.method="FisherChiSquare",enrichFile="$outDir/$id_base.difsea")
}
date()
q('no')

#;

open O,">","$outDir/$id_base.R";
print O  $rscp;
close O;

`R CMD BATCH $outDir/$id_base.R $outDir/$id_base.Rout`;


