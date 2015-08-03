#!/usr/bin/perl -w 
#Author: Pengcheng Yang
#Email: yangpc@mail.biols.ac.cn

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname); 

my $usage= "\n$0 
--gmt    \t<str>\tannotation file in gmt format
--minG    \t<str>\tminimal genes number of one gene sets (default 5)
--maxG    \t<str>\tmaximal genes number of one gene sets (default 500)
--rat     \t<str>\tthe overlap ratio between parent and child, the parent will
\t\tbe removed if greater than this ratio. (default 0.9)
--outDir  \t<str>\toutput directory

NOTE: this script is used to reduce the redundancy of GO terms by their annotation
\tand hierarchy structures
\n";

my($gmt,$outDir,$minG,$maxG,$rat);

GetOptions(
		"gmt:s"=>\$gmt,
		"outDir:s"=>\$outDir,
		"minG:i"=>\$minG,
		"maxG:i"=>\$maxG,
		"rat:f"=>\$rat
);

die $usage if !defined $gmt;

$minG ||= 5;
$maxG ||= 500;
$rat ||= 0.9;

`mkdir -p $outDir` if ! -e $outDir;
my $gmt_base=basename $gmt;

my $rscp = qq#
library(graph)
library(RBGL)
library(AnnotationForge)
library(GO.db)
Childs <- c("GO:0003674","GO:0008150","GO:0005575")

tmp <- scan("$gmt",sep="\\n",what="character")
go2g <- sapply(tmp,function(x){
  x1 <- strsplit(x,"\\t")[[1]]
  x1[3:length(x1)]
})
names(go2g) <- sapply(tmp,function(x) strsplit(x,"\\t")[[1]][1])
go2g <- go2g[setdiff(names(go2g),c("all",Childs))]
go2g <- go2g[intersect(names(go2g),keys(GO.db))]

goids.tb <- select(GO.db,keys=names(go2g),keytypes="GOID",column=c("GOID","ONTOLOGY"))
goids.bp <- goids.tb[goids.tb[,2]=="BP",1]
goids.mf <- goids.tb[goids.tb[,2]=="MF",1]
goids.cc <- goids.tb[goids.tb[,2]=="CC",1]
goids.aces <- c(mget(goids.bp,GOBPANCESTOR),mget(goids.mf,GOMFANCESTOR),mget(goids.cc,GOCCANCESTOR))
goids.off <- c(mget(goids.bp,GOBPOFFSPRING),mget(goids.mf,GOMFOFFSPRING),mget(goids.cc,GOCCOFFSPRING))
goids.child <- c(mget(goids.bp,GOBPCHILDREN),mget(goids.mf,GOMFCHILDREN),mget(goids.cc,GOCCCHILDREN))

goids <- unique(c(unlist(goids.aces),unlist(goids.off),names(go2g)))
goids <- goids[-(grep("all",goids))]
goids <- goids[-(which(is.na(goids)))]
goids <- as.vector(goids)
goterm <- mget(goids,GOTERM)

goLev <- rep("a",length=length(goterm))
names(goLev) <- names(goterm)
lev <- 1
while(length(Childs)>=1){
	loc <- match(Childs,names(goLev))
	goLev[loc][goLev[loc] == "a"] <- lev
	Childs <- unique(as.character(unlist(goids.child[Childs])))
	Childs <- Childs[! is.na(Childs)]
	lev <- lev+1
}

go2g.tb <- as.matrix(data.frame(Gene=as.character(unlist(go2g)),
		GO=rep(names(go2g),sapply(go2g,length))))

go2g.net <- ftM2graphNEL(go2g.tb[,2:1], edgemode="directed")
go.net.bp <- makeGOGraph(ont = "bp")
go.net.mf <- makeGOGraph(ont = "mf")
go.net.cc <- makeGOGraph(ont = "cc")
\#get all the parents from current nodes

bp.ids <- intersect(nodes(go.net.bp),nodes(go2g.net))
mf.ids <- intersect(nodes(go.net.mf),nodes(go2g.net))
cc.ids <- intersect(nodes(go.net.cc),nodes(go2g.net))
go.net.bp.sub <- subGraph(bp.ids,go.net.bp)
go.net.mf.sub <- subGraph(mf.ids,go.net.mf)
go.net.cc.sub <- subGraph(cc.ids,go.net.cc)
save.image(file="$outDir/EACO.RData")

reduce.go <- function(gonet,go2gNet){
	\#calculate the statistics of every term
	nr <- numNodes(gonet)
	go.dt <- data.frame(
			GO = nodes(gonet),
			Ngenes=rep(0,nr),
			GOclass = rep("NA",nr),
			nChild = rep(0,nr),
			nParent = rep(0,nr),
			MchildGO = rep("NA",nr),
			Mchild = rep(0,nr),
			MparentGO = rep("NA",nr),
			Mparent = rep("NA",nr),
			Level=goLev[nodes(gonet)],
			Term= as.character(sapply(nodes(gonet),function(x) Term(goterm[[x]]))),
			stringsAsFactors=FALSE
			)
	
	rownames(go.dt) <- nodes(gonet)
	for(i in nodes(gonet)){
		genes <- edges(go2gNet,i)[[1]]
		go.dt[i,"Ngenes"] <- length(genes)
		go.dt[i,"GOclass"] <- "CC"
		Childs <- inEdges(i,gonet)[[1]]
		Parents <- edges(gonet,i)[[1]]
		go.dt[i,"nChild"] <- length(Childs)
		go.dt[i,"nParent"] <- length(Parents)
		\#for children
		if(length(Childs)>0){
			Child.stat <- sapply(Childs,function(x){
				genesC <- edges(go2gNet,x)[[1]]
				NgenesC <- length(genesC)
				length(intersect(genes,genesC))/go.dt[i,"Ngenes"]
			})
			go.dt[i,"MchildGO"] <- names(Child.stat)[which.max(Child.stat)]
			go.dt[i,"Mchild"] <- max(Child.stat)
		}
		\#for parents
		if(length(Parents)>0){
			Parent.stat <- sapply(Parents,function(x){
					genesP <- edges(go2gNet,x)[[1]]
					NgenesP <- length(genesP)
					length(intersect(genes,genesP))/NgenesP
					})
			go.dt[i,"MparentGO"] <- names(Parent.stat)[which.max(Parent.stat)]
			go.dt[i,"Mparent"] <- max(Parent.stat)
		}
		
	}
	
	filt.stat <- rep(0,6)
	names(filt.stat) <- c("totalSets","Sets_minG","Sets_maxG","distincG","distincG_minG","distincG_maxG")
	filt.stat["totalSets"] <- dim(go.dt)[1]
	filt.stat["distincG"] <- length(unique(as.character(unlist(edges(go2gNet,go.dt\$GO)))))
	go.dt.filt <- go.dt[go.dt\$Ngenes >= $minG,]
	filt.stat["Sets_minG"] <- dim(go.dt.filt)[1]
	filt.stat["distincG_minG"] <- length(unique(as.character(unlist(edges(go2gNet,go.dt.filt\$GO)))))
	go.dt.filt <- go.dt.filt[go.dt.filt\$Ngenes <= $maxG,];dim(go.dt.filt)
	filt.stat["Sets_maxG"] <- dim(go.dt.filt)[1]
	filt.stat["distincG_maxG"] <- length(unique(as.character(unlist(edges(go2gNet,go.dt.filt\$GO)))))
	sub.sub <- subGraph(rownames(go.dt.filt),gonet)
	sub.edge.stat <- data.frame(ChildN=sapply(inEdges(sub.sub),length),
			ParentN=sapply(edges(sub.sub),length))

	\#test ratio gradient
	rt.brk <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1)
	rt.dt <- matrix(nr=length(rt.brk),nc=2)
	colnames(rt.dt) <- c("SetsNum","distingcG")
	rownames(rt.dt) <- paste("r",rt.brk,sep=".")
	for(rt in rt.brk){
		sub.rt <- sub.sub
		sub.dt <- go.dt.filt
		for(i in rownames(go.dt.filt)){
			if(go.dt.filt[i,"Mchild"] > rt){
				if(is.element(go.dt.filt[i,"MchildGO"],rownames(go.dt.filt))){
					sub.rt <- removeNode(i,sub.rt)
				}
			}
		}
		sub.dt <- sub.dt[nodes(sub.rt),]
		rt.dt[paste("r.",rt,sep=""),"SetsNum"] <- dim(sub.dt)[1]
		rt.dt[paste("r.",rt,sep=""),"distingcG"] <- length(unique(as.character(unlist(edges(go2gNet,sub.dt\$GO)))))
	}

	res <- list(dt=go.dt,dt.filt=go.dt.filt,net=sub.sub,stat=filt.stat,ratcut=rt.dt)
	res
}

net.bp.filt <- reduce.go(go.net.bp.sub,go2g.net)
net.mf.filt <- reduce.go(go.net.mf.sub,go2g.net)
net.cc.filt <- reduce.go(go.net.cc.sub,go2g.net)
save(net.bp.filt,net.mf.filt,net.cc.filt,file="$outDir/AfterFilt.RData")
res.dt.filt <- rbind(net.bp.filt[["dt.filt"]],net.mf.filt[["dt.filt"]],net.cc.filt[["dt.filt"]])
res.dt <- rbind(net.bp.filt[["dt"]],net.mf.filt[["dt"]],net.cc.filt[["dt"]])
res.dt.dif <- res.dt[setdiff(rownames(res.dt),rownames(res.dt.filt)),]

\#total filtering statistics
stat.tot <- rep(0,6)
names(stat.tot) <- names(net.bp.filt[["stat"]])
stat.tot["totalSets"] <- dim(res.dt)[1]
stat.tot["distincG"] <- length(unique(as.character(unlist(go2g))))
stat.tot["Sets_minG"] <- dim(res.dt[res.dt[,"Ngenes"]>=$minG,])[1]
stat.tot["Sets_maxG"] <- dim(res.dt[res.dt[,"Ngenes"]>=$minG & res.dt[,"Ngenes"]<$maxG,])[1]
stat.tot["distincG_minG"] <- length(unique(as.character(unlist(go2g[res.dt[res.dt[,"Ngenes"]>=$minG,1]]))))
stat.tot["distincG_maxG"] <- length(unique(as.character(unlist(go2g[res.dt[res.dt[,"Ngenes"]>=$minG & res.dt[,"Ngenes"]<$maxG,1]]))))
res.stat <- rbind(net.bp.filt[["stat"]],net.mf.filt[["stat"]],net.cc.filt[["stat"]],stat.tot)
rownames(res.stat) <- c("BP","MF","CC","total")

rt.brk <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1)
res.rt <- cbind(net.bp.filt[["ratcut"]],net.mf.filt[["ratcut"]],net.cc.filt[["ratcut"]])
colnames(res.rt) <- paste(rep(c("BP","MF","CC"),each=2),colnames(res.rt),sep=".")
res.rt <- data.frame(Ratio=rt.brk,res.rt)
write.table(res.stat,file="$outDir/$gmt_base.filt.stat.txt",sep="\\t",quote=FALSE,col.names=NA)
write.table(res.dt,file="$outDir/$gmt_base.beforefilt.txt",sep="\\t",quote=FALSE,row.names=FALSE)
write.table(res.rt,file="$outDir/$gmt_base.rat.gradient.txt",sep="\\t",quote=FALSE,row.names=FALSE)
q('no')

#;

open O,">","$outDir/$gmt_base.R";
print O $rscp;
close O;

`R CMD BATCH $outDir/$gmt_base.R $outDir/$gmt_base.Rout`;

