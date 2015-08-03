#!/usr/bin/perl -w 
#Author: Pengcheng Yang
#Email: yangpc@mail.biols.ac.cn

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname); 

my $usage= "\n$0 
--gmt     \t<str>\tgene set file in gmt format
--Etb     \t<str>\tenrichment analysis results table. every column is the fdr/pv value of one
\t\tsample, every line is one gene sets, the first column is gene set ID
--Etbcfg  \t<str>\tconfigure file for Etb, in the format:
\t\tSampID<tb>Direction<tb>Groupid...
\t\tSampID: is the columan head of --Etb file
\t\tDirection: is the Up Down of the input
--epv     \t<str>\tpvalue cutoff to visualize (default 0.05)
--outDir  \t<str>\toutput directory
--step    \t<str>\tsteps:
\t\t1. filtering the Etb based on epv and convert the enrichment value into -log10
\t\t2. clustering the gsets based on the enrichment value
\t\t3. visualize the filtered enrichment value based on the order of the clustering
--kapmt   \t<str>\tkappa statistics in matrix format. used to cluster the gSets

\n";

my($gmt,$Etb,$Etbcfg,$epv,$outDir,$step,$kapmt);

GetOptions(
		"gmt:s"=>\$gmt,
		"Etb:s"=>\$Etb,
		"Etbcfg:s"=>\$Etbcfg,
		"epv:f"=>\$epv,
		"outDir:s"=>\$outDir,
		"step:s"=>\$step,
		"kapmt:s"=>\$kapmt
);

die $usage if !defined $Etb;

my $Etb_base = basename $Etb;
my $dir = dirname $0;
$epv ||= 0.05;

`mkdir -p $outDir` if ! -e "$outDir";

my($rscp);

if($step =~ /1/){
	$rscp=qq#
date()
\#filtering the enrichment results
options(stringsAsFactors=FALSE)
et <- read.table("$Etb",header=TRUE,sep="\\t",check.names=FALSE)
rownames(et) <- toupper(et[,1])
et <- as.matrix(et[,-1])
et[is.na(et)] <- 1
et[et==0] <- 1e-5
etcfg <- read.table("$Etbcfg",sep="\\t",check.names=FALSE)
et <- et[,etcfg[,1]]
\#filt enrichment result
et.sig <- et[apply(et,1,function(x) any(x<$epv)),]
et.sig.log <- -log10(et.sig)
write.table(et.sig.log,file="$outDir/$Etb_base.logp",sep="\\t",quote=FALSE)

date()
q('no')
#;

	open O,">","$outDir/$Etb_base.sig.R";
	print O $rscp;
	close O;

	open O,">","$outDir/step1_filt_Etb.sh";
	print O "R CMD BATCH $outDir/$Etb_base.sig.R $outDir/$Etb_base.sig.Rout\n";
	close O;

	`R CMD BATCH $outDir/$Etb_base.sig.R $outDir/$Etb_base.sig.Rout`;

}

if($step =~ /2/){
	open O,">","$outDir/step2_clustering_Etb.sh";
	`perl $dir/run_WGCNA.pl --expr $outDir/$Etb_base.logp --expCut 0 --TOMType signed  --outDir $outDir/01.wgcna --step 1`;
	`perl $dir/run_WGCNA.pl --expr $outDir/$Etb_base.logp --expCut 0 --TOMType signed  --outDir $outDir/01.wgcna --step 2 --beta 8`;
	print O "perl $dir/run_WGCNA.pl --expr $outDir/$Etb_base.logp --expCut 0 --TOMType signed  --outDir $outDir/01.wgcna --step 12 --beta 8 \n\n";
	print O "perl $dir/cluster_based_on_numeric_matrix.pl --mat $outDir/$Etb_base.logp --outDir $outDir/02.cluster/ --TOMdis $outDir/01.wgcna/$Etb_base.logp.networkConstruction.TOM.RData --test\n\n";
	close O;
}

if($step =~ /3/){
	`mkdir -p $outDir/03.vis` if ! -e "$outDir/03.vis";
	$rscp = qq#
options(stringsAsFactors=FALSE)
library(ggplot2)
library(gridExtra)
library(ggdendro)
library(dynamicTreeCut)
library(reshape2)

\#read in data sets
load("$outDir/02.cluster/$Etb_base.logp.tree.RData")
kp <- read.table("$kapmt",check.names=FALSE)
kp <- 1-kp
colnames(kp) <- toupper(colnames(kp))
rownames(kp) <- toupper(rownames(kp))
etcfg <- read.table("$Etbcfg",se="\\t")
et.sig.log <- read.table("$outDir/$Etb_base.logp",sep="\\t")
kp.sig <- kp[rownames(et.sig.log),rownames(et.sig.log)]

tmp <- scan("$gmt",what="character",sep="\\n")
gsets <- sapply(tmp,function(x) strsplit(x,"\\t")[[1]][-c(1,2)])
names(gsets) <- sapply(tmp,function(x) strsplit(x,"\\t")[[1]][1])
setDes <- sapply(tmp,function(x) strsplit(x,"\\t")[[1]][2])
names(setDes) <- toupper(names(gsets))
names(gsets) <- toupper(names(gsets))

\#clustering of gene sets based on kappa statistics
hcl.lst.al[["KP"]] <- hclust(as.dist(kp.sig),method="average")
dynamic.tree[["KP"]] <- cutreeDynamic(dendro=hcl.lst.al[["KP"]],minClusterSize=2,distM=as.matrix(kp.sig),
        cutHeight=as.numeric(quantile(hcl.lst.al[["KP"]]\$height,0.90)))
names(dynamic.tree[["KP"]]) <- hcl.lst.al[["KP"]]\$label

\#output the gene set description of every gene set group
dynamic.tree.dt <- cbind(as.data.frame(dynamic.tree),Desc=setDes[names(dynamic.tree[[1]])])
write.csv(dynamic.tree.dt,file="$outDir/03.vis/$Etb_base.cluster_desc.csv")

\#calculate and output the adjusted rand index between every partition pairs
library(clusterSim)
adjrdi <- function(x,y){
	x1 <- signif(comparing.Partitions(x,y,type="crand"),digit=2) \#rand index
	x2 <- signif(comparing.Partitions(x,y,type="rand"),digit=2) \#adjusted rand index
	text(max(x)*0.1,max(y)*0.7,paste(x1,"\\n",x2,sep=""),pos=4)
}


barp <- function(x){
	a <- barplot(table(x),plot=F)
	rect(as.numeric(a-1),0,as.numeric(a+1),0.9*max(x)*table(x)/max(table(x)))
}

pdf("$outDir/03.vis/$Etb_base.tree_pairs_indx.pdf",width=10,height=10)
pairs(as.data.frame(dynamic.tree),lower.panel=adjrdi,diag.panel=barp,gap=0)
dev.off()


dirct <- rep(etcfg[,2],each=NROW(et.sig.log))
c3 <- rep(etcfg[,3],each=NROW(et.sig.log))
go.dt <- data.frame(term=rep(rownames(et.sig.log),NCOL(et.sig.log)),
		log.p=as.numeric(as.matrix(et.sig.log)),dirct=dirct,c3=c3)
go.dt\$c3 <- factor(go.dt\$c3,levels=sort(unique(etcfg[,3])))
if(dim(etcfg)[2] >3){
	c4 <- rep(etcfg[,4],each=NROW(et.sig.log))
	go.dt <- data.frame(go.dt,c4=c4)
	go.dt\$c4 <- factor(go.dt\$c4,levels=sort(unique(etcfg[,4])))
}
if(dim(etcfg)[2] > 4){
	c5 <- rep(etcfg[,5],each=NROW(et.sig.log))
	go.dt <- data.frame(go.dt,c5=c5)
	go.dt\$c5 <- factor(go.dt\$c5,levels=sort(unique(etcfg[,5])))
}


godt2 <- 0
if(length(grep("Down",etcfg[,2]))>=1){
	drct.sign <- sapply(as.character(go.dt\$dirct),function(x) ifelse(identical(x,"Down"),-1,1))
	dt.nam <- setdiff(names(go.dt),c("log.p","dirct"))
	x.past <- apply(go.dt[,dt.nam],1,function(x) paste(x,collapse="@"))
	logp2 <- tapply(go.dt\$log.p*drct.sign,x.past,function(x){
			x[which.max(abs(x))]
			})
	godt2 <- t(sapply(names(logp2),function(x) strsplit(x,"@")[[1]]))
	colnames(godt2) <- dt.nam
	godt2 <- data.frame(godt2,log.p=logp2)
	rownames(godt2) <- NULL
	et.sig.log2 <- dcast(godt2,term~c3+c4)
	rownames(et.sig.log2) <- et.sig.log2[,1]
	et.sig.log2 <- et.sig.log2[,-1]
	write.table(et.sig.log2,file="$outDir/03.vis/$Etb_base.logp.merg",sep="\\t",quote=FALSE,col.names=NA)
}

setDes.short <- sapply(setDes,function(x) strsplit(x,"\\\\|")[[1]][1])
save(hcl.lst.al,etcfg,dynamic.tree,abr,dst.lst,go.dt,godt2,setDes,setDes.short,et.sig.log,et.sig.log2,file="$outDir/03.vis/$Etb_base.logp.vis.RData")

go.dt[go.dt[,2]>4,2] <- 4

\#plot the top 10 sets
setids <- apply(et.sig.log2,2,function(x) {
		x1 <- which((length(x)-rank(abs(x),ties.method="first"))<11)
		rownames(et.sig.log2)[x1]
	})
setids <- unique(as.vector(setids))
godt.n <- go.dt
godt.n\$term <- factor(go.dt\$term,levels=unique(as.vector(setids)))
godt.n <- godt.n[!is.na(godt.n\$term),]
godt.n.hcl <- hclust(dist(et.sig.log[setids,]),method="average")

write.table(data.frame(
			et.sig.log[godt.n.hcl\$order,],setDes=setDes[setids][godt.n.hcl\$order]),
		file="$outDir/03.vis/$Etb_base.vis.top10.xls",sep="\\t",col.names=NA,quote=FALSE)

p.bar <- ggplot(godt.n,aes(term,log.p,fill=dirct))+
	facet_grid(.~c3)+
	geom_bar(position="dodge",width=0.8,stat="identity")+
	scale_fill_manual("dirct",values=c("red","blue"),limits=c("Up","Down"))+
	theme(axis.text.x = element_text(colour="black"),
			axis.text.y=element_text(colour="black"),
			plot.margin=unit(c(2,2,2,-1),"lines"))+
	ylab("-log10(AdjustedPv)")+xlab("")+
	scale_x_discrete(limits=setids[godt.n.hcl\$order],labels=substr(setDes.short[setids][godt.n.hcl\$order],1,40))
pdf("$outDir/03.vis/$Etb_base.vis.top10.pdf")
print(p.bar + coord_flip(ylim=c(0,4)))
dev.off()


\#plot every cluster tree
for(i in names(hcl.lst.al)){
	et.dhcl <- as.dendrogram(hcl.lst.al[[i]])
	et.ddata <- dendro_data(et.dhcl, type="rectangle")
	dev.off()
	\#setord <- rownames(et.sig.log)[hcl.lst.al[[i]]\$order]
	setord <- names(sort(dynamic.tree[[i]]))
	p.bar <- ggplot(go.dt,aes(term,log.p,fill=dirct))+
		facet_grid(.~c3+c4)+
		geom_bar(position="dodge",width=0.8,stat="identity")+
		scale_fill_manual("dirct",values=c("red","blue"),limits=c("Up","Down"))+
		theme(axis.text.x = element_text(colour="black"),
				axis.text.y=element_text(colour=ifelse(sort(dynamic.tree[[i]]) %% 2 ==1,"black","red")),
				plot.margin=unit(c(2,2,2,-1),"lines"))+
		\#  ylim(0,max(go.dt\$log.p)+1)+
		ylab("-log10(AdjustedPv)")+xlab("")+
		coord_flip(ylim=c(0,4))+
		\#  scale_x_discrete(limits=rev(setord)) +
		scale_x_discrete(limits=setord,labels=substr(setDes.short[setord],1,40))
		
		\#add rectangle for every gene sets
	cl.tb <- c(0,cumsum(table(sort(dynamic.tree[[i]]))))
	for(j in 1:(length(cl.tb)-1)){
		i.col <- ifelse(j %% 2 ==0,"lightblue","antiquewhite")
		if(j %% 2 ==0) next
		rec.dt <- data.frame(xmin=cl.tb[j]+0.5,xmax=cl.tb[j+1]+0.5,ymin=-Inf,ymax=Inf)
		p.bar <- p.bar + geom_rect(data=rec.dt,aes(NULL,NULL,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
				fill=i.col)
	}
	p.bar <- p.bar+
		geom_hline(aes(yintercept=seq(0,4,1)),color="white") +
		geom_vline(aes(xintercept=1:NROW(et.sig.log)),color="white") +
		geom_bar(position="dodge",width=0.8,stat="identity") 

  \#add background color for strip
	p.bar.grb <- ggplotGrob(p.bar)
\#	idx <- 0
\#	cols <- rep("gray",66)
\#	cols[seq(1,66,3)] <- rep(c("red","blue"),each=11)
\#	cols[seq(2,66,3)] <- rep(rep(c("red","blue"),c(6,5)),2)
\#	for( g in 1:length(p.bar.grb\$grobs) ){
\#		if( grepl( "strip.absoluteGrob" , p.bar.grb\$grobs[[g]]\$name ) ){
\#			idx <- idx + 1
\#				sb <- which( grepl( "strip\\\\.background" , names( p.bar.grb\$grobs[[g]]\$children ) ) )
\#				p.bar.grb\$grobs[[g]]\$children[[sb]]\$gp\$fill <- cols[idx]
\#		}
\#	}

  pdf(file=paste("$outDir/03.vis/enrich_logp_gsets_hcl.",i,".pdf",sep=""),
			height=(10+dim(et.sig.log)[1])/6,width=(5+dim(et.sig.log)[2])/2)
	grid.draw(p.bar.grb)
	\#print(p.bar)
	dev.off()

	p.dendro <- ggplot(segment(et.ddata)) +
		geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
		coord_flip() + scale_y_reverse(expand=c(0,0)) + scale_x_continuous(expand=c(0,0.3)) +
		xlab("") + ylab("Height") +
		theme(panel.grid.major.y = element_blank(),
				panel.grid.minor.y = element_blank(),
				axis.text.y = element_blank(),
				axis.ticks.y = element_blank(),
				plot.margin=unit(c(dim(etcfg)[2],-1,2,2),"lines"))
		\#geom_text(data=et.ddata\$labels,aes(x=x, y=y, label=label))
	p.text <-ggplot(data.frame(x=rep(1,length(setDes.short)),y=1:length(setDes.short))) +
		geom_text(aes(x=x,y=y,lab=setDes.short))
}


\#heat map visualization

\#heatmap.2 visualization
library(gplots)
sig2 <- dcast(godt2,term~c3)
rownames(sig2) <- sig2[,1]
sig2 <- sig2[,-1]
sig2 <- as.matrix(sig2)
sig2.log <- log2(abs(sig2)+1)
sig2.log <- sig2.log*sign(sig2)
sig2.log.rnm <- rownames(sig2.log)
rownames(sig2.log) <- as.character(sapply(setDes.short[rownames(sig2.log)],function(x) substr(x,1,40)))
pdf("$outDir/03.vis/$Etb_base.heatmap.2.pdf",width=8,height=25)
hv <- heatmap.2(sig2.log,cexCol=2,cexRow=1,keysize=1,col=bluered, tracecol="\#303030",
		Colv=FALSE,dendrogram="row",
		key.par=list(mar=c(3,4,20,1)),
		key.title="",margins=c(6,20))
dev.off()

out.carpet <- data.frame(t(hv\$carpet),setDes=setDes[sig2.log.rnm[hv\$rowInd]])
rownames(out.carpet) <- sig2.log.rnm[hv\$rowInd]
out.carpet <- out.carpet[dim(out.carpet)[1]:1,]

save(out.carpet,sig2.log,setDes,setDes.short,file="$outDir/03.vis/$Etb_base.heatmap.2.RData")
write.table(out.carpet,file="$outDir/03.vis/$Etb_base.heatmap.2.carpet.xls",sep="\\t",quote=FALSE,col.names=NA)

date()
q('no')
#;

	open O,">","$outDir/03.vis/$Etb_base.vis.R";
	print O $rscp;
	close O;

#`R CMD BATCH $outDir/$Etb_base.vis.R $outDir/$Etb_base.vis.Rout`;

}



