#!/usr/bin/perl -w 
#Author: Pengcheng Yang
#Email: yangpc@mail.biols.ac.cn

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname); 

my $usage= "\n$0 
--etb1    \t<str>\tenrichment result table
--etb2    \t<str>\tthe second enrichment result table to compare to the first etb
--gmt     \t<str>\tgene set data in gmt format
--outDir  \t<str>\toutput direcotry

NOTE: this script is used to plot the Circos graph
\n";

my($etb1,$etb2,$gmt,$outDir);

GetOptions(
		"etb1:s"=>\$etb1,
		"etb2:s"=>\$etb2,
		"gmt:s"=>\$gmt,
		"outDir:s"=>\$outDir
);

die $usage if !defined $etb1;

my $dir=dirname $0;
my $flg_et2="N";
if(!defined $etb2){
	$flg_et2="N";
	$etb2="N";
}else{
	$flg_et2="T";
}

my $rscp=qq#
library(RCircos)
options(stringsAsFactors=FALSE)
source("$dir/EACO_R_Fun.R")
chrlen <- 1000
etb1 <- read.table("$etb1",check.names=FALSE,row.names=1)
dt <- data.frame(Chromosome=rep(rownames(etb1),each=2),
		chromStart=rep(c(1,1+chrlen/2),dim(etb1)[1]),
		chromEnd=rep(c(chrlen/2,chrlen),dim(etb1)[1]))
htg <- t(as.data.frame(apply(etb1,1,function(x) list(t(matrix(x,nr=2))))))
colnames(htg) <- gsub(".Up|.Down","",colnames(etb1)[seq(1,dim(etb1)[2],2)])
etb1.hist <- data.frame(dt,htg)
etb1.cyt <- data.frame(Chromosome=rep(rownames(etb1),each=2),
		ChromStart=rep(c(1,1+chrlen/2),dim(etb1)[1]),
		chromEnd=rep(c(chrlen/2,chrlen),dim(etb1)[1]),
		Band=rep(c("qA1","qC1"),dim(etb1)[1]),
		Stain=rep(c("gneg","gpos"),dim(etb1)[1]))

if(identical("T","$flg_et2")){
etb2 <- read.table("$etb2",check.names=FALSE,row.names=1)
dt <- data.frame(Chromosome=rep(rownames(etb2),each=2),
		chromStart=rep(c(1,1+chrlen/2),dim(etb2)[1]),
		chromEnd=rep(c(chrlen/2,chrlen),dim(etb2)[1]))
htg <- t(as.data.frame(apply(etb2,1,function(x) list(t(matrix(x,nr=2))))))
colnames(htg) <- gsub(".Up|.Down","",colnames(etb2)[seq(1,dim(etb2)[2],2)])
etb2.hist <- data.frame(dt,htg)
etb2.cyt <- data.frame(Chromosome=rep(rownames(etb2),each=2),
		ChromStart=rep(c(1,1+chrlen/2),dim(etb2)[1]),
		chromEnd=rep(c(chrlen/2,chrlen),dim(etb2)[1]),
		Band=rep(c("qA1","qC1"),dim(etb2)[1]),
		Stain=rep(c("gneg","gpos"),dim(etb2)[1]))
}

gmt <- read.gmt("$gmt")

gos.comb <- t(combn(names(gmt[["gsets"]]),2))
gos.comb <- rbind(gos.comb,matrix(rep(names(gmt[["gsets"]]),2),nc=2))
geneN <- length(unique(unlist(gmt[["gsets"]])))
kval <- apply(gos.comb,1,function(x){
		C11 <- length(intersect(gmt[["gsets"]][[x[1]]],gmt[["gsets"]][[x[2]]]))
		C00 <- geneN - length(unique(union(gmt[["gsets"]][[x[1]]],gmt[["gsets"]][[x[2]]])))
		Cg1 <- length(gmt[["gsets"]][[x[1]]])
		Cg2 <- length(gmt[["gsets"]][[x[2]]])
		C10 <- Cg1 - C11
		C01 <- Cg2 - C11
		Oab <- (C11 + C00) / geneN
		Aab <- (Cg1 * Cg2 + (geneN - Cg1) * (geneN - Cg2))/(geneN*geneN)
		Kab <- (Oab - Aab)/ (1- Aab)
		c(Kab,range(c(C11/Cg1,C11/Cg2)))
	})
kval <- round(t(kval),digit=5)

res <- data.frame(gos.comb,kval)
res <- res[res[,3]>0.8,]
res <- rbind(res,res[res[,1] != res[,2],c(2,1,3,4,5)])
colnames(res) <- c("Set1","Set2","Kappa","RatMin","RatMax")



#;

open O,">","$outDir/Rcircos.R";
print O $rscp;
close O;

