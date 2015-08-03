#!/usr/bin/perl -w 
#Author: Pengcheng Yang
#Email: yangpc@mail.biols.ac.cn

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname); 

my $usage= "\n$0 
--etb     \t<str>\tenrichment result table in the format column one sample, row one gene set
--gmt     \t<str>\tgene sets file in gmt format
--outDir  \t<str>\toutput directory

--pv      \t<str>\tpvalue cutoff to filter etb file. at least one sample greater than this value
\t\t(default 0.05)

--ovlpR   \t<str>\toverlap ratio cutoff (default 0.1).
\t\tthe overlap ratio is calculated by overlap_gene_num/gene_num_in_gsets. 
--minOvlp \t<str>\tminimal overlap numbers (default 3)

--mcls    \t<str>\tmethod class. currently SEA, GSEA and MEA (default GSEA)

--rnklist \t<str>\tgene file list, in the format:
\t\tID<tb>FilePath
\t\tfor the file, every file in the format:
\t\tGeneID<tb>Pvalue<tb>logFC

#for SEA based enrichment
--idlist  \t<str>\tgene id file list, in the format:
\t\tSampleID<tb>FilePath
\t\tfor the file, every file contain gene ids with one id one line.

--method  \t<str>\tenrichment analysis method.
\tfor GSEA currently support:
\t\tfrom piano R package: fisher, stouffer, reporter, tailStrength, wilcoxon, mean, median, sum, maxmean, gsea or page.
\t\tevery method can be concatenated by comma. (default is gsea)
\tfor SEA currently support:
\t\tFisherChiSquare, HyperGeometric

--bakid   \t<str>\tfor SEA. background id file, one id one line, if undefined the all
\t\tids in the gmt will be used as background.


NOTE: this program adopt the algorithm from GOFunction to filter
\tthe enrichment results according two criteria:
\tlocal redundancy and global redundancy
\n";

my($etb,$pv,$gmt,$outDir,$method,$rnklist,$mcls,$idlist,$ovlpR,$minOvlp,$bakid);

GetOptions(
		"etb:s"=>\$etb,
		"gmt:s"=>\$gmt,
		"outDir:s"=>\$outDir,
		"method:s"=>\$method,
		"mcls:s"=>\$mcls,
		"rnklist:s"=>\$rnklist,
		"pv:f"=>\$pv,
		"ovlpR:f"=>\$ovlpR,
		"minOvlp:i"=>\$minOvlp,
		"idlist:s"=>\$idlist,
		"bakid:s"=>\$bakid
		);

die $usage if !defined $gmt;

$pv ||= 0.05;
$ovlpR ||= 0.1;
$minOvlp ||= 3;

my $etb_base = basename $etb;
my $gmt_base = basename $gmt;
my $dir = dirname $0;

`mkdir -p $outDir` if ! -e "$outDir";

#produce the local and global gene sets
`mkdir -p $outDir/GO_local` if ! -e "$outDir/GO_local";
`mkdir -p $outDir/Global_local` if ! -e "$outDir/Global_local";
`mkdir -p $outDir/Global` if ! -e "$outDir/Global";
`mkdir -p $outDir/Combin` if ! -e "$outDir/Combin";

	my $rscp = qq#
options(stringsAsFactors=FALSE)
library(GO.db)
\#read in database
\#gene sets
tmp <- scan("$gmt",what="character",sep="\\n")
gsets <- sapply(tmp,function(x) strsplit(x,"\\t")[[1]][-c(1,2)])
names(gsets) <- sapply(tmp,function(x) strsplit(x,"\\t")[[1]][1])
setDes <- sapply(tmp,function(x) strsplit(x,"\\t")[[1]][2])
names(setDes) <- toupper(names(gsets))
names(gsets) <- toupper(names(gsets))

etb <- read.table("$etb",sep="\\t",header=TRUE,row.names=1,check.names=FALSE)
etb <- etb[which(apply(etb,1,function(x) any(x<$pv))),]
rownames(etb) <- toupper(rownames(etb))
write.table(etb,file="$outDir/orig.$etb_base.filtpv",sep="\\t",quote=FALSE,col.names=NA)

gsets <- gsets[rownames(etb)]
setDes <- setDes[rownames(etb)]

goids <- names(gsets)[grep("^GO:",names(gsets))]
goids.onto <- select(GO.db,keys=goids,columns=c("GOID","ONTOLOGY"),keytype="GOID")
goids.all <- tapply(goids.onto[,1],goids.onto[,2],as.character)
goids.off <- list()
goids.off <- c(goids.off,mget(goids.all[["BP"]],GOBPOFFSPRING))
goids.off <- c(goids.off,mget(goids.all[["MF"]],GOMFOFFSPRING))
goids.off <- c(goids.off,mget(goids.all[["CC"]],GOCCOFFSPRING))
goids.off <- sapply(goids.off,function(x) intersect(x,goids))
goids.off <- goids.off[sapply(goids.off,length)>0]

\#local gene sets
off.dif.gmt <- data.frame(setid="a",des="b",gids="c")
off.dif.tb <- data.frame(ancesID="a",offID="b",ancesGnum=1,offGnum=1,spcNum=0)
for(i in names(goids.off)){
	i.in <- intersect(goids.off[[i]],goids)
	if(length(i.in)==0) next
	offgids <- unique(as.character(unlist(gsets[i.in])))
	spcid <- setdiff(gsets[[i]],offgids)
	if(length(spcid)>2)
		off.dif.gmt <- rbind(off.dif.gmt,c(
					i,setDes[i],
					paste(spcid,collapse="\\t")))
	off.dif.tb <- rbind(off.dif.tb,c(
				i,paste(i.in,collapse=","),length(gsets[[i]]),length(offgids),
				length(spcid)))
}
off.dif.gmt <- off.dif.gmt[-1,]
off.dif.tb <- off.dif.tb[-1,]
write.table(off.dif.gmt,file="$outDir/GO_local/local.gmt",sep="\\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(off.dif.tb,file="$outDir/GO_local/local.gmt.gNum",sep="\\t",quote=FALSE,row.names=FALSE)

\#global gene sets
setids <- names(gsets)
glb.gmt.tb <- data.frame(set1="a",set2="b",OverlapN=1,set1_spcN=1,set2_spcN=2,set1_spcI="a",set2_spcI="b")
glb.gmt <- data.frame(setid="a",des="b",gids="a")
for(i in 1:(length(gsets)-1)){
	for(j in (i+1):length(gsets)){
		if(is.element(setids[j],goids.off[[setids[i]]]) || 
				is.element(setids[i],goids.off[[setids[j]]])) next
		commid <- intersect(gsets[[setids[i]]],gsets[[setids[j]]])
		if(length(commid) < $minOvlp) next
		if(length(commid)/length(gsets[[setids[j]]]) < $ovlpR &&
				length(commid)/length(gsets[[setids[j]]]) < $ovlpR) next
		spcid.i <- setdiff(gsets[[setids[i]]],gsets[[setids[j]]])
		spcid.j <- setdiff(gsets[[setids[j]]],gsets[[setids[i]]])
		glb.gmt.tb <- rbind(glb.gmt.tb,
				c(setids[i],setids[j],
					length(commid),
					length(spcid.i),
					length(spcid.j),
					paste(spcid.i,collapse=","),
					paste(spcid.j,collapse=",")))
		if(length(spcid.i)>2){
			glb.gmt <- rbind(glb.gmt,
					c(paste(setids[i],"_SPC_VS_",setids[j],sep=""),
						paste(setids[i],"_SPC_VS_",setids[j],sep=""),
						paste(spcid.i,collapse="\\t")))
		}

		if(length(spcid.j)>2){
			glb.gmt <- rbind(glb.gmt,
					c(paste(setids[i],"_VS_",setids[j],"_SPC",sep=""),
						paste(setids[i],"_VS_",setids[j],"_SPC",sep=""),
						paste(spcid.j,collapse="\\t")))
		}

	}
}
glb.gmt <- glb.gmt[-1,]
glb.gmt.tb <- glb.gmt.tb[-1,] 
write.table(glb.gmt,file="$outDir/Global/glb.gmt",sep="\\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(glb.gmt.tb,file="$outDir/Global/glb.gmt.tb",sep="\\t",quote=FALSE,row.names=FALSE)

\#global local sets
ids.off <- vector("list",length=length(gsets))
names(ids.off) <- names(gsets)
ids.off[names(goids.off)] <- goids.off
for(i in 1:dim(glb.gmt.tb)[1]){
	if(glb.gmt.tb[i,4]>2 && glb.gmt.tb[i,5]<3){
		ids.off[[glb.gmt.tb[i,1]]] <- append(ids.off[[glb.gmt.tb[i,1]]],glb.gmt.tb[i,2])
	}
	if(glb.gmt.tb[i,5]>2 && glb.gmt.tb[i,4]<3){
		ids.off[[glb.gmt.tb[i,2]]] <- append(ids.off[[glb.gmt.tb[i,2]]],glb.gmt.tb[i,1])
	}
}

glocal.gmt <- data.frame(setid="a",des="a",gids="a")
glocal.tb <- data.frame(ancesID="a",offID="a",acesGnum=1,offGnum=1,spcNum=1)
ids.off <- ids.off[sapply(ids.off,length)>0]
ids.off <- sapply(ids.off,function(x) unique(x))
for(i in names(ids.off)){
	offgids <- unique(as.character(unlist(gsets[ids.off[[i]]])))
	spcid <- setdiff(gsets[[i]],offgids)
	if(length(spcid)>2)
		glocal.gmt <- rbind(glocal.gmt,c(
					paste(i,"_glocal",sep=""),setDes[i],
					paste(spcid,collapse="\\t")))
	glocal.tb <- rbind(glocal.tb,c(
				paste(i,"_glocal",sep=""),paste(ids.off[[i]],collapse=","),length(gsets[[i]]),length(offgids),
				length(spcid)))
}
glocal.gmt <- glocal.gmt[-1,]
glocal.tb <- glocal.tb[-1,]
write.table(glocal.gmt,file="$outDir/Global_local/glocal.gmt",sep="\\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(glocal.tb,file="$outDir/Global_local/glocal.tb",sep="\\t",quote=FALSE,row.names=FALSE)


write.table(rbind(off.dif.gmt,glb.gmt,glocal.gmt),file="$outDir/Global/glb.local.gmt",sep="\\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

q('no')
#;
open O,">","$outDir/GOFunction.produceGMT.GSEA.R";
print O $rscp;
close O;

if(! -e "$outDir/Global/glb.local.gmt"){
	`R CMD BATCH $outDir/GOFunction.produceGMT.GSEA.R $outDir/GOFunction.produceGMT.GSEA.Rout`;
}

if($mcls eq "GSEA"){

#run piano do enrichment
  open O,">","$outDir/01.run_piano.sh";
  open O2,">","$outDir/Etb.cfg";
  open I,$rnklist;
	my(@info,$id);
  while(<I>){
    chomp;
    @info=split /\t/;
    my $rnk_base = basename $info[1];
    print O "perl $dir/run_piano.pl --gmt $outDir/Global/glb.local.gmt --outDir $outDir/$rnk_base.enrich --glevel $info[1] --method $method --minG 3 --maxG 1950\n";
    print O2 "$rnk_base.Up\tDown\n$rnk_base.Down\tUp\n"; #the second column Up/Down denote the direction in B for the format A-VS-B
  }
  close I; close O; close O2;

#parse piano results
  open O2,">","$outDir/02.parse_piano.sh";
	foreach $id (split /,/,$method){
		open I1,$rnklist;
		open O,">","$outDir/enrich.$id";
		while(<I1>){
			chomp;@info=split /\t/;
			my $rnk_base = basename $info[1];
			print O "$info[0]\t$outDir/$rnk_base.enrich/$info[0].piano.$id\n";
		} 
		close I1; close O;
		print O2 "\n#parse piano file\n";
		print O2 "perl $dir/parse_piano.pl --outDir $outDir/ --pianof $outDir/enrich.$id\n\n";
	}
  close O2;

#local and global redunction
	$rscp = qq#
\#read in data sets
options(stringsAsFactors=FALSE)
source("$dir/EACO_R_Fun.R")

orig.enrich <- read.table("$outDir/orig.enrich.$method.padj.tb.filtpv",sep="\\t",header=TRUE,row.names=1,check.names=FALSE)
sub.enrich <- read.table("$outDir/enrich.$method.pv.tb",sep="\\t",header=TRUE,row.names=1,check.names=FALSE)
sub.enrich[is.na(sub.enrich)] <- 1

go.local.gmt <- read.gmt("$outDir/GO_local/local.gmt")
glb.gmt <- read.gmt("$outDir/Global/glb.gmt")
glocal.gmt <- read.gmt("$outDir/Global_local/glocal.gmt")
orig.gmt <- read.gmt("$gmt")

go.local.tb <- read.table("$outDir/GO_local/local.gmt.gNum",sep="\\t",header=TRUE,row.names=1,check.names=FALSE)
glb.tb <- read.table("$outDir/Global/glb.gmt.tb",sep="\\t",header=TRUE,check.names=FALSE)
glocal.tb <- read.table("$outDir/Global_local/glocal.tb",sep="\\t",header=TRUE,row.names=1,check.names=FALSE)
rownames(glocal.tb) <- gsub("_glocal","",rownames(glocal.tb))

\#GO local redundancy
local.res <- GOFunction.local(go.local.tb,orig.enrich,sub.enrich,orig.gmt,pv=$pv,out.dir="$outDir/GO_local")
write.table(orig.enrich[names(local.res[["setDes"]]),],file="$outDir/GO_local/$etb_base.local.tb",sep="\\t",col.names=NA,
		quote=FALSE)

\#global redundancy
glb.res <- GOFunction.glb(glb.tb,orig.enrich,sub.enrich,orig.gmt,pv=$pv,out.dir="$outDir/Global")
write.table(orig.enrich[names(glb.res[["setDes"]]),],file="$outDir/Global/$etb_base.global.tb",sep="\\t",col.names=NA,
		quote=FALSE)

\#Global local redundancy
sub.enrich.glocal <- sub.enrich[grep("glocal",rownames(sub.enrich)),]
rownames(sub.enrich.glocal) <- gsub("_glocal","",rownames(sub.enrich.glocal))
glocal.res <- GOFunction.local(glocal.tb,orig.enrich,sub.enrich.glocal,orig.gmt,pv=$pv,out.dir="$outDir/Global_local")
write.table(orig.enrich[names(glocal.res[["setDes"]]),],file="$outDir/Global_local/$etb_base.glocal.tb",sep="\\t",col.names=NA,
		quote=FALSE)

\#combine the results together and output
\#intersect of these three results
commids <- intersect(intersect(names(glb.res[["setDes"]]),names(local.res[["setDes"]])),names(glocal.res[["setDes"]]))
mergids <- unique(c(names(glb.res[["setDes"]]),names(local.res[["setDes"]]),names(glocal.res[["setDes"]])))
merg.des <- data.frame(gsetID=mergids,
		go.local=rep(NA,length(mergids)),
		global=rep(NA,length(mergids)),
		glocal=rep(NA,length(mergids)))
rownames(merg.des) <- mergids
merg.des[names(local.res[["setDes"]]),"go.local"] <- local.res[["setDes"]]
merg.des[names(glb.res[["setDes"]]),"global"] <- glb.res[["setDes"]]
merg.des[names(glocal.res[["setDes"]]),"glocal"] <- glocal.res[["setDes"]]
write(commids,file="$outDir/Combin/gset.commids")
write(mergids,file="$outDir/Combin/gset.mergids")
write.table(merg.des,file="$outDir/Combin/combin.setDes.tb",row.names=FALSE,quote=FALSE,sep="\\t")

write.table(orig.enrich[mergids,],file="$outDir/Combin/$etb_base.mergids.tb",sep="\\t",col.names=NA,quote=FALSE)
write.table(orig.enrich[commids,],file="$outDir/Combin/$etb_base.commids.tb",sep="\\t",col.names=NA,quote=FALSE)

save(local.res,glb.res,glocal.res,file="$outDir/GOFunction_reduction.RData")
q('no')

#;

	#check whether the piano was run
	if(-e "$outDir/enrich.$method.pv.tb"){
		open O,">","$outDir/GOFunction.reduce.GSEA.R";
		print O $rscp;
		close O;
		if(! -e "$outDir/GOFunction_reduction.RData"){
			`R CMD BATCH $outDir/GOFunction.reduce.GSEA.R $outDir/GOFunction.reduce.GSEA.Rout`;
		}
	}

}

if($mcls eq "SEA"){
	
#run SEA enrichment for the preconstructed gene sets
	my(@info);
	open I,$idlist;
	open O,">","$outDir/01.run_SEA.sh";
	open O2,">","$outDir/difsea.list";
	while(<I>){
		chomp;
		@info=split /\t/;
		my $file_base = basename $info[1];
		print O "perl $dir/run_SEA.pl --gmt $outDir/Global/glb.local.gmt --outDir $outDir/01.Enrich --supplyID $info[1] --bakid $bakid\n";
		print O2 "$outDir/01.Enrich/$file_base.difsea\n";
	}
	close O;close O2;

#cat the enrichment together
	open O,">","$outDir/02.pase_SEA.sh";
	print O "#cat together the enrichment results\n";
	print O "perl $dir/cat_enrich_result_to_pvtb.pl --difcfg $outDir/difsea.list --outDir $outDir/\n\n";
	close O;

#local and global redunction
	$rscp = qq#
\#read in data sets
options(stringsAsFactors=FALSE)
source("$dir/EACO_R_Fun.R")

orig.enrich <- read.table("$outDir/orig.difsea.list.qv.tb.filtpv",sep="\\t",header=TRUE,row.names=1,check.names=FALSE)
sub.enrich <- read.table("$outDir/difsea.list.pv.tb",sep="\\t",header=TRUE,row.names=1,check.names=FALSE)

go.local.gmt <- read.gmt("$outDir/GO_local/local.gmt")
glb.gmt <- read.gmt("$outDir/Global/glb.gmt")
glocal.gmt <- read.gmt("$outDir/Global_local/glocal.gmt")
orig.gmt <- read.gmt("$gmt")

go.local.tb <- read.table("$outDir/GO_local/local.gmt.gNum",sep="\\t",header=TRUE,row.names=1,check.names=FALSE)
glb.tb <- read.table("$outDir/Global/glb.gmt.tb",sep="\\t",header=TRUE,check.names=FALSE)
glocal.tb <- read.table("$outDir/Global_local/glocal.tb",sep="\\t",header=TRUE,row.names=1,check.names=FALSE)
rownames(glocal.tb) <- gsub("_glocal","",rownames(glocal.tb))

\#GO local redundancy
local.res <- GOFunction.local(go.local.tb,orig.enrich,sub.enrich,orig.gmt,pv=$pv,out.dir="$outDir/GO_local")
write.table(orig.enrich[names(local.res[["setDes"]]),],file="$outDir/GO_local/$etb_base.local.tb",sep="\\t",col.names=NA,
		quote=FALSE)

\#global redundancy
glb.res <- GOFunction.glb(glb.tb,orig.enrich,sub.enrich,orig.gmt,pv=$pv,out.dir="$outDir/Global")
write.table(orig.enrich[names(glb.res[["setDes"]]),],file="$outDir/Global/$etb_base.global.tb",sep="\\t",col.names=NA,
		quote=FALSE)

\#Global local redundancy
sub.enrich.glocal <- sub.enrich[grep("glocal",rownames(sub.enrich)),]
rownames(sub.enrich.glocal) <- gsub("_glocal","",rownames(sub.enrich.glocal))
glocal.res <- GOFunction.local(glocal.tb,orig.enrich,sub.enrich.glocal,orig.gmt,pv=$pv,out.dir="$outDir/Global_local")
write.table(orig.enrich[names(glocal.res[["setDes"]]),],file="$outDir/Global_local/$etb_base.glocal.tb",sep="\\t",col.names=NA,
		quote=FALSE)

\#combine the results together and output
\#intersect of these three results
commids <- intersect(intersect(names(glb.res[["setDes"]]),names(local.res[["setDes"]])),names(glocal.res[["setDes"]]))
mergids <- unique(c(names(glb.res[["setDes"]]),names(local.res[["setDes"]]),names(glocal.res[["setDes"]])))
merg.des <- data.frame(gsetID=mergids,
		go.local=rep(NA,length(mergids)),
		global=rep(NA,length(mergids)),
		glocal=rep(NA,length(mergids)))
rownames(merg.des) <- mergids
merg.des[names(local.res[["setDes"]]),"go.local"] <- local.res[["setDes"]]
merg.des[names(glb.res[["setDes"]]),"global"] <- glb.res[["setDes"]]
merg.des[names(glocal.res[["setDes"]]),"glocal"] <- glocal.res[["setDes"]]
write(commids,file="$outDir/Combin/gset.commids")
write(mergids,file="$outDir/Combin/gset.mergids")
write.table(merg.des,file="$outDir/Combin/combin.setDes.tb",row.names=FALSE,quote=FALSE,sep="\\t")

write.table(orig.enrich[mergids,],file="$outDir/Combin/$etb_base.mergids.tb",sep="\\t",col.names=NA,quote=FALSE)
write.table(orig.enrich[commids,],file="$outDir/Combin/$etb_base.commids.tb",sep="\\t",col.names=NA,quote=FALSE)

save(local.res,glb.res,glocal.res,file="$outDir/GOFunction_reduction.RData")
q('no')

#;

	#check whether the piano was run
	if(-e "$outDir/difsea.list.pv.tb"){
		open O,">","$outDir/GOFunction.reduce.GSEA.R";
		print O $rscp;
		close O;
		if(! -e "$outDir/GOFunction_reduction.RData"){
			`R CMD BATCH $outDir/GOFunction.reduce.GSEA.R $outDir/GOFunction.reduce.GSEA.Rout`;
		}
	}

}





