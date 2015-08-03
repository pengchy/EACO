#!/usr/bin/perl -w 
#Author: Pengcheng Yang
#Email: yangpc@mail.biols.ac.cn

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname); 

my $usage= "\n$0 
--inf     \t<str>\tinput file
--informt \t<str>\tinput file format:
\t\twego: geneID<tb>GOID<tb>GOID<tb>...
\t\tthree: geneID<tb>ClassID<tb>ClassDesc
\t\tgaf: GO annotation file format
--outDir  \t<str>\toutput directory
\n";

my($inf,$informt,$outDir);

GetOptions(
		"inf:s"=>\$inf,
		"informt:s"=>\$informt,
		"outDir:s"=>\$outDir
);

die $usage if !defined $informt;

my $inf_base = basename $inf;

my(@info,%gen2set,%set2gen,%setDes,$rscp);
if($informt eq "wego"){
	$rscp = qq#
date()
library(GO.db)
goterm <- as.list(GOTERM)
tmp <- scan("$inf",what="character",sep="\\n")

g2go <- sapply(tmp,function(x){
		x1 <- strsplit(x,"\\t")[[1]]
		x1[2:length(x1)]
	})
names(g2go) <- sapply(tmp,function(x) strsplit(x,"\\t")[[1]][1])

go2g <- tapply(rep(names(g2go),sapply(g2go,length)),unlist(g2go),as.character)
go2g <- go2g[intersect(names(go2g),names(goterm))]

gmt <- sapply(names(go2g),function(x) paste(x,Term(goterm[[x]]),paste(go2g[[x]],collapse="\\t"),sep="\\t"))
write.table(gmt,file="$outDir/$inf_base.gmt",quote=FALSE,row.names=FALSE,col.names=FALSE)

date()
q('no')
#;

	open O,">","$outDir/$inf_base.R";
	print O $rscp;
	close O;

`R CMD BATCH $outDir/$inf_base.R $outDir/$inf_base.Rout`;

}

if($informt eq "three"){
	$rscp=qq#
date()
options(stringsAsFactors=FALSE)
tmp <- read.table("$inf",sep="\\t",check.names=FALSE,quote="")
tmp <- tmp[grep("IPR",tmp[,2]),]
set2g <- tapply(tmp[,1],tmp[,2],as.character)
setDes <- unique(paste(tmp[,2],tmp[,3],sep="\\t"))
names(setDes) <- sapply(setDes,function(x) strsplit(x,"\\t")[[1]][1])
setDes <- sapply(setDes,function(x) gsub("IPR[0-9]+\\t","",x))
gmt <- sapply(names(set2g),function(x) paste(x,setDes[x],paste(set2g[[x]],collapse="\\t"),sep="\\t"))
write.table(gmt,file="$outDir/$inf_base.gmt",col.names=FALSE,row.names=FALSE,quote=FALSE)
date()
q('no')
#;

	open O,">","$outDir/$inf_base.R";
	print O $rscp;
	close O;

	`R CMD BATCH $outDir/$inf_base.R $outDir/$inf_base.Rout`;

}

if($informt eq "gaf"){
	my(%set2gene,%seqDes);
	open I,$inf;
	while(<I>){
		chomp;@info=split /\t/;
		next if /^!/;
		$set2gene{$info[3]}{$info[1]}++;
		$seqDes{$info[3]}++;
	}
	close I;

	open O,">","$outDir/$inf_base.gmt";
	close O;
}


