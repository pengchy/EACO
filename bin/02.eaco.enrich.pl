#!/usr/bin/perl -w 
#Author: Pengcheng Yang
#Email: yangpc@mail.biols.ac.cn

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname); 

my $usage= "\n$0 
--gmt     \t<str>\tgene set file in the gmt format
--outDir  \t<str>\toutput directory
--minG    \t<str>\tallowed minimal gene number of every gene sets (default 5)
--maxG    \t<str>\tallowed maximal gene number of every gene sets (default 500)
--mcls    \t<str>\tmethod class. currently SEA, GSEA and MEA (default GSEA)
--method  \t<str>\tenrichment analysis method.
\tfor GSEA currently support:
\t\tfrom piano R package: fisher, stouffer, reporter, tailStrength, wilcoxon, mean, median, sum, maxmean, gsea or page.
\t\tevery method can be concatenated by comma. (default is gsea)
\tfor SEA currently support:
\t\tFisherChiSquare, HyperGeometric

#for rank based enrichment
--rnklist \t<str>\tgene file list, in the format:
\t\tID<tb>FilePath
\t\tfor the file, every file in the format:
\t\tGeneID<tb>Pvalue<tb>logFC

#for SEA based enrichment
--idlist  \t<str>\tgene id file list, in the format:
\t\tSampleID<tb>FilePath
\t\tfor the file, every file contain gene ids with one id one line.

--bakid   \t<str>\tfor SEA. background id file, one id one line, if undefined the all
\t\tids in the gmt will be used as background.

NOTE: this script is used to run enrichment analysis using the gmt as gene sets input
\n";

my($gmt,$outDir,$method,$rnklist,$mcls,$minG,$maxG,$idlist,$bakid);

GetOptions(
		"gmt:s"=>\$gmt,
		"outDir:s"=>\$outDir,
		"method:s"=>\$method,
		"mcls:s"=>\$mcls,
		"rnklist:s"=>\$rnklist,
		"minG:i"=>\$minG,
		"maxG:i"=>\$maxG,
		"idlist:s"=>\$idlist,
		"bakid:s"=>\$bakid
);

die $usage if !defined $gmt;

$mcls ||= "GSEA";
$method ||= "gsea";
$minG ||= 5;
$maxG ||= 500;

my $dir = dirname $0;
$outDir = "$outDir/$mcls/";

my(@info,$id);

`mkdir -p $outDir` if ! -e "$outDir";

if($mcls eq "GSEA"){

	open O,">","$outDir/01.run_piano.sh";
	open O2,">","$outDir/Etb.cfg";
	open I,$rnklist;
	while(<I>){
		chomp;
		@info=split /\t/;
		my $rnk_base = basename $info[1];
		print O "perl $dir/run_piano.pl --gmt $gmt --outDir $outDir/$rnk_base.enrich --glevel $info[1] --method $method --minG $minG --maxG $maxG\n";
		print O2 "$rnk_base.Up\tDown\n$rnk_base.Down\tUp\n"; #the second column Up/Down denote the direction in B for the format A-VS-B
	}
	close I; close O; close O2;

	open O2,">","$outDir/02.parse_piano.sh";
	foreach $id (split /,/,$method){
		open I1,$rnklist;
		open O,">","$outDir/enrich.$id";
		while(<I1>){
			chomp;@info=split /\t/;
			my $rnk_base = basename $info[1];
			print O "$info[0]\t$outDir/$rnk_base.enrich/$info[0].piano.$id\n";
		}
		close I1; 
		close O;
		print O2 "\n#parse piano file\n";
		print O2 "perl $dir/parse_piano.pl --outDir $outDir/ --pianof $outDir/enrich.$id\n\n";
		print O2 "perl $dir/run_GOFunction.pl --mcls $mcls --gmt $gmt --outDir $outDir/GOFunction_$id --etb $outDir/enrich.$id.padj.tb --rnklist $rnklist --method $id \n\n";
	}
	close O2; 

}

if($mcls eq "SEA"){
	open I,$idlist;
	open O,">","$outDir/01.run_SEA.sh";
	open O2,">","$outDir/difsea.list";
	open O3,">","$outDir/Etb.cfg";
	`mkdir -p $outDir/01.Enrich` if ! -e "$outDir/01.Enrich";
	while(<I>){
		chomp;
		@info=split /\t/;
		my $file_base = basename $info[1];
		print O "perl $dir/run_SEA.pl --gmt $gmt --supplyID $info[1] --outDir $outDir/01.Enrich --method $method\n";
		print O2 "$outDir/01.Enrich/$file_base.difsea\n";
		print O3 "$file_base.difsea\n";
	}
	close O;close I;close O2;close O3;

	open O,">","$outDir/02.parse_SEA.sh";
	print O "#cat together the enrichment results\n";
	print O "perl $dir/cat_enrich_result_to_pvtb.pl --difcfg $outDir/difsea.list --outDir $outDir/\n\n";
	print O "#run GOfunction to reduce the redundancy\n";
	print O "perl $dir/run_GOFunction.pl --mcls $mcls  --gmt $gmt  --outDir $outDir/GOFunction  --etb $outDir/difsea.list.qv.tb --idlist $idlist --bakid $bakid\n\n";
	close O;
}


