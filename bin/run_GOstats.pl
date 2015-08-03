#!/usr/bin/perl -w 
#Author: Pengcheng Yang
#Email: yangpc@mail.biols.ac.cn

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname); 

my $usage= qq#\n$0 
--gmt     \t<str>\tgene sets file in gmt format
--outDir  \t<str>\toutput directory
--supplyID\t<str>\tsupplied ID to do enrichment. one id one line
--bakid   \t<str>\tbackground id file, one id one line, if undefined the all
\t\tids in the gmt will be used as background.
--padjM   \t<str>p adjust method: "holm","hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"
\t\t(default fdr)

--cls     \t<str>one of GO|KEGG

\n#;

my($gmt,$outDir,$supplyID,$bakid,$padjM,$cls);

GetOptions(
		"gmt:s"=>\$gmt,
		"outDir:s"=>\$outDir,
		"supplyID:s"=>\$supplyID,
		"bakid:s"=>\$bakid,
		"padjM:s"=>\$padjM,
		"cls:s"=>\$cls
);

die $usage if !defined $gmt && !defined $supplyID;

$bakid ||= "N";

my $rscp=qq#
date()
library(GOstats) 
library("GSEABase")
options(stringsAsFactors=FALSE)
tmp <- scan("$gmt",what="character",sep="\\n")
gsets <- sapply(tmp,function(x) strsplit(x,"\\t")[[1]][-c(1,2)])
names(gsets) <- sapply(tmp,function(x) strsplit(x,"\\t")[[1]][1])
setDes <- sapply(tmp,function(x) strsplit(x,"\\t")[[1]][2])
names(setDes) <- names(gsets)

if(identical("$bakid","N")){
	universe <- unique(as.character(unlist(gsets)))
}else{
	universe <- scan("$bakid",sep="\\n",what="character")
}

suppid <- scan("$supplyID",sep="\\n",what="character")

if(identical("$cls","GO")){
	goframeData <- data.frame(frame.go_id=rep(names(gsets),sapply(gsets,length)),
			frame.Evidence=rep("IEA",length(unlist(gsets))),
			frame.gene_id=as.character(unlist(gsets)))
	goFrame=GOFrame(goframeData,organism="Locusta migratoria")
	goAllFrame=GOAllFrame(goFrame)

	gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
	params <- GSEAGOHyperGParams(name="Custom GSEA based annot Params",
			geneSetCollection=gsc,
			geneIds = suppid,
			universeGeneIds = universe,
			ontology = "MF",
			pvalueCutoff = 0.05,
			conditional = FALSE,
			testDirection = "over")
	Over <- hyperGTest(params)

}

#;

open O,">","$outDir/GOstats.R";
print O "$rscp\n";
close O;




