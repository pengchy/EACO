#!/usr/bin/perl -w 
#Author: Pengcheng Yang
#Email: yangpc@mail.biols.ac.cn

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname); 

my $usage= "\n$0 
--genedesc\t<str>\tgene description file, in the format: geneid<tb>description
--gids    \t<str>\tall gene ids of the organism, one gene id one line
--go      \t<str>\tgene ontology annotition file in gmt format
--path    \t<str>\tpathway gene sets in gmt format
--domain  \t<str>\tdomain gene sets in gmt format
--spc     \t<str>\tTwo letters abbreviation of species
--System  \t<str>\tsystem code like Ensembl
--SysAbr  \t<str>\tSystem abbreviation
--outDir  \t<str>\toutput directory
--step    \t<str>\t
\t\t1. primary gene ids
\t\t2. GO relations
\t\t3. KEGG relations

NOTE: this script is used to produce files for GO-Elite program for construct
\torganism specific database

\n";

my($genedesc,$gids,$go,$path,$domain,$outDir,$step,$spc,$System,$SysAbr);

GetOptions(
		"genedesc:s"=>\$genedesc,
		"gids:s"=>\$gids,
		"domain:s"=>\$domain,
		"outDir:s"=>\$outDir,
		"go:s"=>\$go,
		"path:s"=>\$path,
		"spc:s"=>\$spc,
		"System:s"=>\$System,
		"SysAbr:s"=>\$SysAbr,
		"step:s"=>\$step
);

die $usage if !defined $genedesc;

$step ||= "1";


$outDir = "$outDir/EnsMart00/$spc/";
`mkdir -p $outDir` if ! -e "$outDir";

if($step =~ /1/){
	`mkdir $outDir/gene` if ! -e "$outDir/gene";
	`perl -e 'print "ID\\tSymbol\\tDescription\\n"' > $outDir/gene/$System.txt`;
	`cat $genedesc |perl -ne 'chomp;\@info=split /\\t/;print "\$info[0]\\t\$info[0]\\t\$info[1]\n"' >> $outDir/gene/$System.txt`;
	print qq#join -t \$'\\t' -v 2 <(sort $genedesc) <(sort $gids) |perl -ne 'chomp;print "\$_\\t\$_\\tNA\\n"' >> $outDir/gene/$System.txt\n#;
}

if($step =~ /2/){
	`mkdir $outDir/gene-go` if ! -e "$outDir/gene-go";
	`mkdir $outDir/nested` if ! -e "$outDir/nested";
	`perl -e 'print "Gene ID\\tGO ID\n"' > $outDir/gene-go/$System-GeneOntology.txt`;
	`cat $go |perl -ne 'chomp;\@info=split /\\t/;map{print "\$_\\t\$info[0]\\n"}\@info[2..\$#info]' >> $outDir/gene-go/$System-GeneOntology.txt`;
	`perl -e 'print "$System\\tontology_id\n"' > $outDir/nested/$System\_to_Nested-GO.txt`;
	`cat $go |perl -ne 'chomp;\@info=split /\\t/;map{print "\$_\\t\$info[0]\\n"}\@info[2..\$#info]' >> $outDir/nested/$System\_to_Nested-GO.txt`;
	
}

if($step =~ /3/){
	`mkdir $outDir/gene-mapp` if ! -e "$outDir/gene-mapp";
	if(defined $path){
		`perl -e 'print "GeneID\\t$System\\tPathway\n"' > $outDir/gene-mapp/$System-KEGG.txt`;
		`cat $path |perl -ne 'chomp;\@info=split /\\t/;map{print "\$_\\t$SysAbr\\t\$info[1]:\$info[0]\\n"}\@info[2..\$#info]' >> $outDir/gene-mapp/$System-KEGG.txt` if defined $path;
	}
	if(defined $domain){
		`perl -e 'print "GeneID\\t$System\\tPathway\n"' > $outDir/gene-mapp/$System-DOMAIN.txt`;
		`cat $domain |perl -ne 'chomp;\@info=split /\\t/;map{print "\$_\\t$SysAbr\\t\$info[1]:\$info[0]\\n"}\@info[2..\$#info]' >> $outDir/gene-mapp/$System-DOMAIN.txt` if defined $domain;
	}
}



