#!/usr/bin/perl -w 
#Author: Pengcheng Yang
#Email: yangpc@mail.biols.ac.cn

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname); 

my $usage= "\n$0 
--gmt     \t<str>\tgmt format gene sets file
--outDir  \t<str>\toutput directory

NOTE: this script is used to plot the correlation between the maxG cut off and 
\tthe distinct gene number. this information will be used to define the maxG cut off
\tfor the reduce_GO_by_anno_hierarch.pl program

the output graph with xaxis the gene number cutoff, yaxis is the percentage of gene set
number with gene number greater than this cutoff, or the percentage of distinct gene number
of the gene set with the gene number less than this cutoff.

\n";

my($gmt,$outDir);

GetOptions(
		"gmt:s"=>\$gmt,
		"outDir:s"=>\$outDir
);

die $usage if !defined $gmt;

my $gmt_base=basename $gmt;

my $rscp=qq#
options(stringsAsFactors=FALSE)
tmp <- scan("$gmt",sep="\\n",what="character")
set2g <- sapply(tmp,function(x){
		x1 <- strsplit(x,"\\t")[[1]]
		x1[3:length(x1)]
		})
names(set2g) <- sapply(tmp,function(x) strsplit(x,"\\t")[[1]][1])
set.name <- sapply(tmp,function(x){
		strsplit(x,"\\t")[[1]][2]})
names(set.name) <- names(set2g)

set2g.len <- sapply(set2g,length)
write.table(data.frame(gSets=names(set2g),
			genNum=as.numeric(set2g.len),
			setNam=as.character(set.name)),
		file="$outDir/$gmt_base.gNum.tb",
		row.names=FALSE,quote=FALSE,sep="\\t")
brks <- seq(5,max(set2g.len)+5,5)
dt <- matrix(0,nr=length(brks),nc=5)
rownames(dt) <- brks
colnames(dt) <- c("brks","SetNum","SetNumPerc","GeneNum","GeneNumPerc")
dt[,"brks"] <- brks
dt[,"SetNum"] <- as.numeric(sapply(brks,function(x){
			length(set2g.len[set2g.len<x])
		}))
dt[,"SetNumPerc"] <- 100*dt[,"SetNum"]/length(set2g.len)
dt[,"GeneNum"] <- as.numeric(sapply(brks,function(x){
			length(unique(as.character(unlist(set2g[set2g.len<x]))))
			}))
dt[,"GeneNumPerc"] <- 100*dt[,"GeneNum"]/length(unique(as.character(unlist(set2g))))

write.table(dt,file="$outDir/$gmt_base.maxG_and_distinctGeneNum.tb",sep="\\t",row.names=FALSE,quote=FALSE)
dt.long <- data.frame(brks=rep(dt[,"brks"],2),
		Percentage=as.numeric(dt[,c("SetNumPerc","GeneNumPerc")]),
		Group=rep(c("SetNumPerc","GeneNumPerc"),each=NROW(dt))
		)

save(dt.long,file="$outDir/$gmt_base.data.RData")
x.brks <- c(1,5,10,20,50,100,200,500,800,1000,1500,2000,2500,3000,6000,10000)
library(ggplot2)
p.line <- ggplot(dt.long,aes(brks,Percentage,colour=Group))+geom_line()+
	scale_x_log10(breaks=x.brks[1:(length(which(x.brks<max(brks)))+1)])+
	xlab("GeneNumber cutoff")+
	theme(axis.text.x = element_text(angle=60))

pdf("$outDir/$gmt_base.maxG_and_distinctGeneNum.tb.plot.pdf",width=10,height=7)
plot(p.line)
dev.off()

q('no')
#;

open O,">","$outDir/maxG_and_DistinctGeneNum.R";
print O $rscp;
close O;

`R CMD BATCH $outDir/maxG_and_DistinctGeneNum.R $outDir/maxG_and_DistinctGeneNum.Rout`;

