#!/usr/bin/perl -w 
#Author: Pengcheng Yang
#Email: yangpc@mail.biols.ac.cn

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname); 

my $usage= "\n$0 
--GO      \t<str>\tdo filtering for GO annotation
--gSets   \t<str>\tgene set file, in the gmt format:
\t\tclass<tb>desc<tb>geneid1<tb>geneid2<tb>...
--outDir  \t<str>\toutput directory

NOTE: this script is used to prepare the gene sets data for enrichment.
For GO annotation data, the filtering will be carried according to the three criteria
minG, maxG and rat of reduce_GO_by_anno_hierarch.pl program

the other no hierarchy gene sets will be filtered based on the three criteria:
minG, maxG and kappa statistics of every two gene sets

all these parameters should be set on the specific directory

\n";

my($GO,$gSets,$outDir);

GetOptions(
		"GO:s"=>\$GO,
		"gSets:s"=>\$gSets,
		"outDir:s"=>\$outDir
);

die $usage if !defined $gSets;

my $dir = dirname $0;
my $gSets_base = basename $gSets;

if(defined $GO){
	`mkdir -p $outDir/GO` if ! -e "$outDir/GO";
	open O,">","$outDir/GO/filtering_GO.log";
	print O "#map the GO ontology to all their ancestors\n";
	print O "perl $dir/map_go_ancestor.pl --gmt $gSets --outDir $outDir/GO/\n\n";
	print O "#test the influence of maxG cutoff on distinct gene number\n";
	print O "perl $dir/maxG_and_DistinctGeneNum.pl --gmt $outDir/GO/$gSets_base.addaces.gmt --outDir $outDir/GO/\n\n";
	print O "perl $dir/reduce_GO_by_anno_hierarch.pl --gmt $outDir/GO/$gSets_base.addaces.gmt --outDir $outDir/GO\n\n";
	print O "#filter the result according to the minG, maxG and ratio defined above\n";
	print O "cat $outDir/GO/$gSets_base.addaces.gmt.beforefilt.txt |awk '\$2>=\$minG && \$2<= \$maxG && \$7<= \$ratio' > $outDir/GO/$gSets_base.addaces.gmt.afterfilt.txt\n\n";
	print O "cat $outDir/GO/$gSets_base.addaces.gmt.afterfilt.txt |cut -f1 |sort |join -t \$'\\t'  -  <(sort $outDir/GO/$gSets_base.addaces.gmt) > $outDir/GO/$gSets_base.addaces.gmt.filt.gmt\n\n";
	print O qq#echo "after filt, the gset number is:" > $outDir/GO/$gSets_base.addaces.gmt.filt.gmt.stat\n#;
	print O qq#wc -l $outDir/GO/$gSets_base.addaces.gmt.filt.gmt >> $gSets_base.addaces.gmt.filt.gmt.stat\n#;
	print O qq#echo "after filt, the distinct gene number is:" >> $gSets_base.addaces.gmt.filt.gmt.stat\n#;
	print O qq#cat $outDir/GO/$gSets_base.addaces.gmt.filt.gmt |perl -ne 'chomp;\@info=split /\\t/;print join("\\n",\@info[2..\$\#info]),"\\n"'|sort -u |wc -l >> $gSets_base.addaces.gmt.filt.gmt.stat\n\n#;
	close O;
}else{
	`mkdir -p $outDir/Kappa/` if ! -e "$outDir/Kappa";
	open O,">","$outDir/Kappa/run_Kappa.log";
	print O "#test the influence of maxG cutoff on distinct gene number\n";
	print O "perl $dir/maxG_and_DistinctGeneNum.pl --gmt $gSets --outDir $outDir/Kappa/\n\n";
	print O "cat $outDir/Kappa/$gSets_base.gNum.tb |awk '\$2>= && \$2<=' |cut -f1 |join -t \$'\\t' <(sort -) <(sort $gSets) > $outDir/Kappa/$gSets_base.filt.gmt\n\n";
	print O qq#echo "after minG and maxG filt, the gset number is:" > $outDir/Kappa/$gSets_base.filt.gmt.stat\n#;
	print O qq#wc -l $outDir/Kappa/$gSets_base.filt.gmt >> $outDir/Kappa/$gSets_base.filt.gmt.stat\n#;
	print O qq#echo "the distinct gene number is:" >> $outDir/Kappa/$gSets_base.filt.gmt.stat\n#;
	print O qq#cat $outDir/Kappa/$gSets_base.filt.gmt |perl -ne 'chomp;\@info=split /\\t/;print join("\\n",\@info[2..\$\#info]),"\\n"'|sort -u |wc -l >> $outDir/Kappa/$gSets_base.filt.gmt.stat\n\n#;

	print O "#calculate kappa statistics for all gene sets pairs\n";
	print O "perl $dir/calculate.cohen.kappa.pl --gmt $outDir/Kappa/$gSets_base.filt.gmt --outDir $outDir/Kappa/\n\n";
	print O "perl $dir/three_column_to_matrix.pl --dt $outDir/Kappa/$gSets_base.filt.gmt.kappa.tb --mt $outDir/Kappa/$gSets_base.kappa.mt\n\n";
	print O "#filtering gene sets based on the kappa statistics\n";
	print O "perl $dir/filt_gsets_based_on_kappa.pl --gmt $outDir/Kappa/$gSets_base.filt.gmt --kap $outDir/Kappa/$gSets_base.filt.gmt.kappa.tb --outDir $outDir/Kappa/\n\n";
	close O;
}



