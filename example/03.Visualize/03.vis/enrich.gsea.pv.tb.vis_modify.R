
options(stringsAsFactors=FALSE)
library(ggplot2)
library(gridExtra)
library(ggdendro)
library(dynamicTreeCut)

#read in data sets
load("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/02.cluster/enrich.gsea.pv.tb.logp.tree.RData")
kp <- read.table("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/01.prep/Kappa/gSets.gmt.filt.gmt.kappa.mt",check.names=FALSE)
kp <- 1-kp
colnames(kp) <- toupper(colnames(kp))
rownames(kp) <- toupper(rownames(kp))
etcfg <- read.table("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/02.Enrich/GSEA/Etb.cfg",se="\t")
et.sig.log <- read.table("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/enrich.gsea.pv.tb.logp",sep="\t")
kp.sig <- kp[rownames(et.sig.log),rownames(et.sig.log)]

tmp <- scan("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/01.prep/Kappa/gSets.gmt.filt.gmt.filtkappa.gmt",what="character",sep="\n")
gsets <- sapply(tmp,function(x) strsplit(x,"\t")[[1]][-c(1,2)])
names(gsets) <- sapply(tmp,function(x) strsplit(x,"\t")[[1]][1])
setDes <- sapply(tmp,function(x) strsplit(x,"\t")[[1]][2])
names(setDes) <- names(gsets)
names(setDes) <- toupper(names(gsets))

#clustering of gene sets based on kappa statistics
hcl.lst.al[["KP"]] <- hclust(as.dist(kp.sig),method="average")
dynamic.tree[["KP"]] <- cutreeDynamic(dendro=hcl.lst.al[["KP"]],minClusterSize=2,distM=as.matrix(kp.sig),
        cutHeight=as.numeric(quantile(hcl.lst.al[["KP"]]$height,0.90)))
names(dynamic.tree[["KP"]]) <- hcl.lst.al[["KP"]]$label

#output the gene set description of every gene set group
dynamic.tree.dt <- cbind(as.data.frame(dynamic.tree),Desc=setDes[names(dynamic.tree[[1]])])
write.csv(dynamic.tree.dt,file="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/03.vis/enrich.gsea.pv.tb.cluster_desc.csv")

#calculate and output the adjusted rand index between every partition pairs
library(clusterSim)
adjrdi <- function(x,y){
	x1 <- signif(comparing.Partitions(x,y,type="crand"),digit=2) #rand index
	x2 <- signif(comparing.Partitions(x,y,type="rand"),digit=2) #adjusted rand index
	text(max(x)*0.1,max(y)*0.7,paste(x1,"\n",x2,sep=""),pos=4)
}


barp <- function(x){
	a <- barplot(table(x),plot=F)
	rect(as.numeric(a-1),0,as.numeric(a+1),0.9*max(x)*table(x)/max(table(x)))
}

pdf("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/03.vis/enrich.gsea.pv.tb.tree_pairs_indx.pdf",width=10,height=10)
pairs(as.data.frame(dynamic.tree),lower.panel=adjrdi,diag.panel=barp,gap=0)
dev.off()

time <- rep(etcfg[,3],each=NROW(et.sig.log))
dirct <- rep(etcfg[,2],each=NROW(et.sig.log))

go.dt <- data.frame(go.term=rep(rownames(et.sig.log),NCOL(et.sig.log)),
		log.p=as.numeric(as.matrix(et.sig.log)),
		time=time,dirct=dirct)
go.dt$time <- factor(go.dt$time,levels=c("egg","1_2","3","4","5","adult"))

save(hcl.lst.al,dynamic.tree,abr,dst.lst,go.dt,setDes,file="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/03.vis/enrich.gsea.pv.tb.logp.vis.RData")

#plot every cluster tree
for(i in names(hcl.lst.al)){
	et.dhcl <- as.dendrogram(hcl.lst.al[[i]])
	et.ddata <- dendro_data(et.dhcl, type="rectangle")
	dev.off()
	#setord <- rownames(et.sig.log)[hcl.lst.al[[i]]$order]
	setord <- names(sort(dynamic.tree[[i]]))
	p.bar <- ggplot(go.dt,aes(go.term,log.p,fill=dirct))+
		facet_grid(.~time)+
		geom_bar(position="dodge",width=0.8,stat="identity")+
		scale_fill_manual("dirct",values=c("red","blue"),limits=c("Up","Down"))+
		theme(axis.text.x = element_text(colour="black"),
				axis.text.y=element_text(colour=ifelse(sort(dynamic.tree[[i]]) %% 2 ==1,"black","red")),
				plot.margin=unit(c(2,2,2,-1),"lines"))+
		#  ylim(0,max(go.dt$log.p)+1)+
		ylim(0,3)+
		ylab("-log10(AdjustedPv)")+xlab("")+
		coord_flip()+
		#  scale_x_discrete(limits=rev(setord)) +
		scale_x_discrete(limits=setord,labels=substr(setDes[setord],1,40))
		
		#add rectangle for every gene sets
	cl.tb <- c(0,cumsum(table(sort(dynamic.tree[[i]]))))
	for(j in 1:(length(cl.tb)-1)){
		i.col <- ifelse(j %% 2 ==0,"lightblue","antiquewhite")
		if(j %% 2 ==0) next
		rec.dt <- data.frame(xmin=cl.tb[j]+0.5,xmax=cl.tb[j+1]+0.5,ymin=-Inf,ymax=Inf)
		p.bar <- p.bar + geom_rect(data=rec.dt,aes(NULL,NULL,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
				fill=i.col)
	}
	p.bar <- p.bar+
		geom_vline(aes(xintercept=1:NROW(et.sig.log)),color="white") +
		geom_hline(aes(yintercept=seq(0,3,1)),color="white") +
		geom_bar(position="dodge",width=0.8,stat="identity") 

  #add background color for strip
	p.bar.grb <- ggplotGrob(p.bar)

  pdf(file=paste("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/03.vis/enrich_logp_gsets_hcl.",i,".pdf",sep=""))
	grid.draw(p.bar.grb)
	#print(p.bar)
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
		#geom_text(data=et.ddata$labels,aes(x=x, y=y, label=label))
	p.text <-ggplot(data.frame(x=rep(1,length(setDes)),y=1:length(setDes))) +
		geom_text(aes(x=x,y=y,lab=setDes))
}

date()
q('no')
