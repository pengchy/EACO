#!/usr/bin/perl -w 
#Author: Pengcheng Yang
#Email: yangpc@mail.biols.ac.cn

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname); 

my $usage= "\n$0 
--gmt     \t<str>\tgo annotation file in gmt format
--outDir  \t<str>\toutput directory

NOTE: this script is used to map all the ancestor of input go terms
\n";

my($gmt,$outDir);

GetOptions(
		"gmt:s"=>\$gmt,
		"outDir:s"=>\$outDir
);

die $usage if !defined $gmt;

my $gmt_base = basename $gmt;

my $rscp = qq#

library(GOstats)

tmp <- scan("$gmt",sep="\\n",what="character")
go2g <- sapply(tmp,function(x){
  x1 <- strsplit(x,"\\t")[[1]]
  x1[3:length(x1)]
})
names(go2g) <- sapply(tmp,function(x) strsplit(x,"\\t")[[1]][1])
go2g <- go2g[intersect(names(go2g),keys(GO.db))]

goframeData <- data.frame(frame.go_id=rep(names(go2g),sapply(go2g,length)),
		frame.Evidence=rep("IEA",length(unlist(go2g))),frame.gene_id=as.character(unlist(go2g)))
goFrame <- GOFrame(goframeData,organism="customized")
goAllFrame <- GOAllFrame(goFrame)
go.data <- getGOFrameData(goAllFrame)
go2g <- tapply(go.data[,3],go.data[,1],as.character)
go2g <- go2g[setdiff(names(go2g),c("all","GO:0003674","GO:0008150","GO:0005575"))]

\#output the gmt file
gmt <- data.frame(GO=names(go2g),
		Term=select(GO.db,keys=names(go2g),keytypes="GOID",column="TERM")[,2],
		genes=as.character(sapply(go2g,function(x) paste(x,collapse="\\t"))))
write.table(gmt,file="$outDir/$gmt_base.addaces.gmt",col.names=FALSE,row.names=FALSE,sep="\\t",quote=FALSE)

res <- c(length(go2g),length(unique(as.character(unlist(go2g)))))
names(res) <- c("totalGenSets","totalGenNum")
write.table(res,file="$outDir/$gmt_base.addaces.gmt.stat",sep="\\t",quote=FALSE)

q('no')
#;

open O,">","$outDir/$gmt_base.addaces.R";
print O $rscp;
close O;

`R CMD BATCH $outDir/$gmt_base.addaces.R $outDir/$gmt_base.addaces.Rout`;



