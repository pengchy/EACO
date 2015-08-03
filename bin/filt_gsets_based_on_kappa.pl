#!/usr/bin/perl -w 
#Author: Pengcheng Yang
#Email: yangpc@mail.biols.ac.cn

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname); 

my $usage= "\n$0 
--gmt     \t<str>\tgene set file in gmt format
--kap     \t<str>\tfile for kappa statistics value of every gset pair in the format
\t\tGset1<tb>Gset2<tb>kappa_value<tb>ratio1<tb>ratio2
--plot    \t<str>\twhether to plot the pairs scatterplot of the kappa, minRat and maxRat of
\t\tevery gene sets pairs. (default NO)
--k       \t<str>\tkappa cutoff, (default 0.97)
--outDir  \t<str>\toutput directory

\n";

my($gmt,$outDir,$kap,$plot,$k);

GetOptions(
		"gmt:s"=>\$gmt,
		"kap:s"=>\$kap,
		"outDir:s"=>\$outDir,
		"plot"=>\$plot,
		"k:f"=>\$k
);

die $usage if !defined $kap;

$k ||= 0.97;

if(defined $plot){
	$plot="YES";
}else{
	$plot="NO";
}

my $kap_base = basename $kap if defined $kap;
my $gmt_base = basename $gmt;

my $rscp = qq#
date()
options(stringsAsFactors=FALSE)
tmp <- scan("$gmt",what="character",sep="\\n")
gsets <- sapply(tmp,function(x) strsplit(x,"\\t")[[1]][-c(1,2)])
names(gsets) <- sapply(tmp,function(x) strsplit(x,"\\t")[[1]][1])
setDes <- sapply(tmp,function(x) strsplit(x,"\\t")[[1]][2])
names(setDes) <- names(gsets)

kap <- read.table("$kap",sep="\\t")
kap <- kap[kap[,1]!=kap[,2],]
colnames(kap) <- c("Gset1","Gset2","Kappa","minRat","maxRat")
if(identical("YES","$plot")){
	pdf("$outDir/$kap_base.kap_min_maxrat.pairs.pdf")
	pairs(kap[,3:5],pch=".")
	dev.off()
}
kap.brk <- seq(0.5,1,0.025)
kap.dt <- matrix(nr=length(kap.brk),nc=3)
colnames(kap.dt) <- c("KappCut","setNum","distincG")
rownames(kap.dt) <- paste("k.",kap.brk,sep="")
kap.dt[,"KappCut"] <- kap.brk
for(k in kap.brk){
	filt.id <- apply(kap[kap[,3]>=k,1:2],1,function(x){
			x1 <- as.character(x)
			x1[which.min(c(length(gsets[[x1[1]]]),length(gsets[[x1[2]]])))]
			})
	remain.id <- setdiff(names(gsets),filt.id)
	kap.dt[paste("k.",k,sep=""),"setNum"] <- length(remain.id)
	kap.dt[paste("k.",k,sep=""),"distincG"] <- length(unique(as.character(unlist(gsets[remain.id]))))
}
write.table(kap.dt,file="$outDir/$kap_base.kappa_gradient.tb",sep="\\t",row.names=FALSE)

\#filtering based on kappa value
library(graph)
dt.filt <- kap[kap[,3]>$k,1:2]
dt.net <- ftM2graphNEL(as.matrix(dt.filt[,1:2]))
net.conn <- connComp(dt.net)
res.gmt <- data.frame(setid="a",desc="a",gids="a")
out.res <- data.frame()
delid <- vector("character")
for(i in 1:length(net.conn)){
	uni.nod <- unique(as.character(unlist(gsets[net.conn[[i]]])))
	tmp.dt <- data.frame(cluster=rep(i,length(net.conn[[i]])),
			ID=net.conn[[i]],
			gNum=as.numeric(sapply(gsets[net.conn[[i]]],length)),
			gNumRatio=signif(as.numeric(sapply(gsets[net.conn[[i]]],length))/length(uni.nod),digits=4),
			Descrp=setDes[net.conn[[i]]])
	rownames(tmp.dt) <- net.conn[[i]]
	if(dim(out.res)[1]==0){
		out.res <- tmp.dt
	}else{
		out.res <- rbind(out.res,tmp.dt)
	}
	slctid <- as.character(tmp.dt[which.max(sapply(tmp.dt[,"Descrp"],nchar)),"ID"])
	uslctid <- setdiff(tmp.dt\$ID,slctid)
	delid <- c(delid,uslctid)
	res.gmt <- rbind(res.gmt,c(
				slctid,
				paste(tmp.dt[slctid,"Descrp"],"|",paste(apply(tmp.dt[uslctid,c("ID","Descrp")],1,function(x){
							x1=as.character(x);paste(x1,collapse=":")}),collapse="|"),sep=""),
				paste(gsets[[slctid]],collapse="\\t")))
}

write.table(out.res,file="$outDir/$gmt_base.kapp.clst",sep="\\t",quote=FALSE,row.names=FALSE)

conn.ids <- as.character(unlist(net.conn))
remain.ids <- setdiff(names(gsets),conn.ids)
res.gmt <- res.gmt[-1,]
res.gmt <- rbind(res.gmt,
		data.frame(setid=remain.ids,
			desc=setDes[remain.ids],
			gids=sapply(gsets[remain.ids],function(x) paste(x,collapse="\\t"))))
write.table(res.gmt,file="$outDir/$gmt_base.filtkappa.gmt",sep="\\t",quote=FALSE,col.names=FALSE,row.names=FALSE)

res.stat <- c(0,0)
names(res.stat) <- c("gSetsNum","gNum")
res.stat["gSetsNum"] <- length(remain.id)
res.stat["gNum"] <- length(unique(as.character(unlist(gsets[remain.id]))))
write.table(res.stat,file="$outDir/$gmt_base.filtkappa.gmt.stat",sep="\\t",quote=FALSE,col.names=FALSE)

date()
q('no')
#;

open O,">","$outDir/$kap_base.filt.R";
print O $rscp;
close O;

`R CMD BATCH $outDir/$kap_base.filt.R $outDir/$kap_base.filt.Rout`;





__END__


mt <- read.table("$mt",sep="\\t",stringsAsFactors=FALSE,check.name=FALSE)
write.table(table(cut(as.numeric(as.matrix(mt)),seq(0,1,0.1),include.lowest=TRUE)),
		file="$outDir/kappa_value_distribute.tb",sep="\\t",quote=FALSE,row.names=FALSE)
mt <- 1-mt
colnames(mt) <- rownames(mt)
mt.dt <- as.dist(mt)
dt.clt <- hclust(mt.dt,method="average")

\#test the effect of height on gene sets
if(! file.exists("$outDir/kappa_gradient")){
	dir.create("$outDir/kappa_gradient")
}

h.brk <- c(0.001,0.01,0.05,0.1,0.15,0.2,0.3)
h.dt <- matrix(nr=length(h.brk),nc=3)
colnames(h.dt) <- c("breaks","setNum","distincG")
rownames(h.dt) <- paste("h.",h.brk,sep="")
h.dt[,"breaks"] <- h.brk
for(h in h.brk){
	clt.cut <- cutree(dt.clt,h=h)
	clt.cut.tb <- table(clt.cut) 
	multi.clt <- names(clt.cut.tb)[which(clt.cut.tb>1)]
	multi.clt <- clt.cut[which(clt.cut %in% multi.clt)]
	multi.clt.dt <- data.frame(gSets=names(multi.clt),cluster=as.numeric(multi.clt),
			gNum=sapply(gsets[names(multi.clt)],length),setDes=setDes[names(multi.clt)])
	write.table(multi.clt.dt,file=paste("$outDir/kappa_gradient/$mt_base.cluster_kappa.tb.",h,".xls",sep=""),sep="\\t",quote=FALSE,row.names=FALSE)
	multi.clt.merg <- t(sapply(tapply(names(multi.clt),as.character(multi.clt),as.character),function(x) {
			x1 <- x[which.max(sapply(gsets[x],length))]
			x2 <- paste(paste(names(setDes[x]),setDes[x],sep="|"),collapse=";")
			c(x1,x2)
			}))
	filt.id <- setdiff(names(multi.clt),multi.clt.merg[,1])
	filt.id.loc <- rownames(mt) %in% filt.id
	remain.id <- rownames(mt)[!filt.id.loc]
	write(remain.id,file=paste("$outDir/kappa_gradient/$mt_base.geset.id.",h,sep=""))
	h.dt[paste("h.",h,sep=""),"setNum"] <- length(filt.id)
	h.dt[paste("h.",h,sep=""),"distincG"] <- length(unique(as.character(unlist(gsets[remain.id]))))
}

write.table(h.dt,file="$outDir/kappa_gradient.stat.tb",sep="\\t",quote=FALSE,row.names=FALSE)
save(h.dt,dt.clt,mt,gsets,setDes,file="$outDir/kappa_gradient.RData")
date()
q('no')
#;


