#!/usr/bin/perl -w 
##Author: Pengcheng Yang
##Email: yangpc@mail.biols.ac.cn
#
use strict;
use Getopt::Long;
use File::Basename qw(basename dirname); 
use GACP qw(parse_config);

my $usage= "\n$0 
--bincfg  \t<str>\tbin configure file, in the format: Name = file_path
--list    \t<str>\tsemsim value list file, in the format: ID<tb>file_path
\t\tthe file_path specify the semsim value file, this file is in the format:
\t\tGO1<tb>GO2<tb>value
--outDir  \t<str>\toutput directory

NOTE: this script is used to plot correlation heatmap of sematic similarity
\n";

my($bincfg,$outDir,$list);

GetOptions(
		"bincfg:s"=>\$bincfg,
		"list:s"=>\$list,
		"outDir:s"=>\$outDir
);

die $usage if !defined $list;

my $list_base=basename $list;
my $outpdf="$outDir/$list_base.cor.pdf";
my $rscp=qq#
tmp <- read.table("$list",stringsAsFactors=FALSE)
dt <- list()
uninam <- vector()
for(i in 1:dim(tmp)[1]){
	tmp1 <- read.table(tmp[i,2],stringsAsFactors=FALSE)
	dt[[tmp[i,1]]] <- as.numeric(tmp1[,3][grep("None",tmp1[,3],invert=TRUE)])
	tmp2 <- paste(tmp1[,1],tmp1[,2],sep="-")
	names(dt[[tmp[i,1]]]) <- tmp2[grep("None",tmp1[,3],invert=TRUE)]
	tmp2 <- tmp2[grep("None",tmp1[,3],invert=TRUE)]
	if(length(uninam)==0)
		uninam <- tmp2
	else
		uninam <- intersect(uninam,tmp2)
}

res <- matrix(0,nr=length(uninam),nc=length(dt))
colnames(res) <- names(dt)
rownames(res) <- uninam
for(i in names(dt)){
	res[,i] <- dt[[i]][uninam]
}

library(gplots)
expr.cor <- cor(res)
pdf("$outpdf",width=8,height=8)
heatmap.2(expr.cor,symm=TRUE,Colv=TRUE,revC=TRUE,trace="none",margins=c(9,9),
dendrogram="both",density.info="none",col=rev(heat.colors(20)),
cellnote=signif(expr.cor,digit=1)
)
dev.off()
#;

open O,">","$outDir/$list_base.R";
print O $rscp;
close O;

`R CMD BATCH $outDir/$list_base.R $outDir/$list_base.Rout`;


