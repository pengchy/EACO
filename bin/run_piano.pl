#!/usr/bin/perl -w 
#Author: Pengcheng Yang
#Email: yangpc@mail.biols.ac.cn

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname); 

my $usage= qq#\n$0 
--glevel  \t<str>\tgene id file to do enrichment, in the format:
\t\tGeneID<tb>Pvalue<tb>logFC<tb>SignedLogPv<tb>rank
\t\tif the file is the format: A_VS_B, the logFC is logFC(B/A)
--minG    \t<str>\tallowed minimal gene number of every gene sets (default 5)
--maxG    \t<str>\tallowed maximal gene number of every gene sets (default 500)
--gmt    \t<str>\tgene set file in gmt format
--method  \t<str>\talgorithm used to calculate the enrichment score, currently support:
\t\tfisher, stouffer, reporter, tailStrength, wilcoxon, mean, median, sum, maxmean, gsea or page.
\t\tevery method can be concatenated by comma

--outDir  \t<str>\toutput directory
NOTE: this script is used to run piano
\n#;

my($outDir,$glevel,$gmt,$method,$minG,$maxG);

GetOptions(
		"glevel:s"=>\$glevel,
		"gmt:s"=>\$gmt,
		"outDir:s"=>\$outDir,
		"method:s"=>\$method,
		"minG:i"=>\$minG,
		"maxG:i"=>\$maxG
);

die $usage if !defined $glevel;
$minG ||= 5;
$maxG ||= 500;

`mkdir $outDir` if ! -e $outDir;
my $glevel_base = basename($glevel);

my $rscp = qq#
date()
library(piano)
\#read in the gene level statistics
tmp = read.table("$glevel",stringsAsFactors=FALSE)
glevel = tmp[,2]
names(glevel) = tmp[,1]
direct = tmp[,4]
names(direct) = tmp[,1]

\#read in gene sets data
tmp = scan("$gmt",what="character",sep="\\n")
tmp2 = sapply(tmp,function(x){
		x1 = strsplit(x,"\\t")[[1]]
		x2 = paste(x1[-c(1,2)],x1[1],sep="\\t")
		x2
		})
tmp2 = unlist(tmp2)
gene2gensets = data.frame(gene=sapply(tmp2,function(x) strsplit(x,"\\t")[[1]][1]),
		set=sapply(tmp2,function(x) strsplit(x,"\\t")[[1]][2]))
row.names(gene2gensets) <- 1:length(tmp2)
gscs = loadGSC(gene2gensets)
if(length(grep("fisher","$method"))==1){
	gsaRes.fisher = runGSA(geneLevelStats=glevel,directions=direct,geneSetStat="fisher",
			signifMethod="geneSampling",adjMethod="fdr",gsc=gscs,gsSizeLim=c($minG,$maxG),ncpus=4)
		GSAsummaryTable(gsaRes.fisher,save=TRUE,file="$outDir/$glevel_base.piano.fisher") 
}

if(length(grep("stouffer","$method"))==1){
	gsaRes.stouffer = runGSA(geneLevelStats=glevel,directions=direct,geneSetStat="stouffer",
			signifMethod="geneSampling",adjMethod="fdr",gsc=gscs,gsSizeLim=c($minG,$maxG),ncpus=4)
		GSAsummaryTable(gsaRes.stouffer,save=TRUE,file="$outDir/$glevel_base.piano.stouffer")
}
 
if(length(grep("reporter","$method"))==1){
	gsaRes.reporter = runGSA(geneLevelStats=glevel,directions=direct,geneSetStat="reporter",
			signifMethod="geneSampling",adjMethod="fdr",gsc=gscs,gsSizeLim=c($minG,$maxG),ncpus=4)
		GSAsummaryTable(gsaRes.reporter,save=TRUE,file="$outDir/$glevel_base.piano.reporter")
}

if(length(grep("tailStrength","$method"))==1){
	gsaRes.tailStrength = runGSA(geneLevelStats=glevel,directions=direct,geneSetStat="tailStrength",
			signifMethod="geneSampling",adjMethod="fdr",gsc=gscs,gsSizeLim=c($minG,$maxG),ncpus=4)
		GSAsummaryTable(gsaRes.tailStrength,save=TRUE,file="$outDir/$glevel_base.piano.tailStrength")
}

if(length(grep("wilcoxon","$method"))==1){
	gsaRes.wilcoxon = runGSA(geneLevelStats=glevel,directions=direct,geneSetStat="wilcoxon",
			signifMethod="geneSampling",adjMethod="fdr",gsc=gscs,gsSizeLim=c($minG,$maxG),ncpus=4)
		GSAsummaryTable(gsaRes.wilcoxon,save=TRUE,file="$outDir/$glevel_base.piano.wilcoxon")
}

if(length(grep("page","$method"))==1){
	gsaRes.page = runGSA(geneLevelStats=direct,geneSetStat="page",
			signifMethod="geneSampling",adjMethod="fdr",gsc=gscs,gsSizeLim=c($minG,$maxG),ncpus=4)
		GSAsummaryTable(gsaRes.page,save=TRUE,file="$outDir/$glevel_base.piano.page")
}

if(length(grep("gsea","$method"))==1){
	gsaRes.gsea = runGSA(geneLevelStats=direct,geneSetStat="gsea",
			signifMethod="geneSampling",adjMethod="fdr",gsc=gscs,gsSizeLim=c($minG,$maxG),ncpus=4)
		GSAsummaryTable(gsaRes.gsea,save=TRUE,file="$outDir/$glevel_base.piano.gsea")
}

if(length(grep("maxmean","$method"))==1){
	gsaRes.maxmean = runGSA(geneLevelStats=direct,geneSetStat="maxmean",
			signifMethod="geneSampling",adjMethod="fdr",gsc=gscs,gsSizeLim=c($minG,$maxG),ncpus=4)
		GSAsummaryTable(gsaRes.maxmean,save=TRUE,file="$outDir/$glevel_base.piano.maxmean")
}

\#save(gsaRes.fisher,gsaRes.gsea,gsaRes.maxmean,gsaRes.page,gsaRes.reporter,gsaRes.stouffer,gsaRes.tailStrength,file="$outDir/$glevel_base.RData")

\#consensus
\#resList <- list(gsaRes.fisher,gsaRes.gsea,gsaRes.maxmean,gsaRes.page,gsaRes.reporter,gsaRes.stouffer,gsaRes.tailStrength)
\#names(resList) <- c("fisher","gsea","maxmean","page","reporter","stouffer","tailStrength")
\#ch <- consensusHeatmap(resList,cutoff=30,method="mean")
date()
q('no')
#;

open O,">","$outDir/piano.R";
print O $rscp;
close O;

`R CMD BATCH $outDir/piano.R $outDir/piano.Rout`;



