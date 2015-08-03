
options(stringsAsFactors=FALSE)
dt <- read.table("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/enrich.gsea.pv.tb.logp",stringsAsFactors=FALSE)
dt <- as.matrix(dt)

#pca representative to judge the cluster numbers
dt.prcomp <- prcomp(t(dt))
dt.pc.load <- as.data.frame(dt.prcomp$rotation)

library(GGally)
library(ggplot2)
aa <- ggpairs(dt.pc.load[,1:3],
              lower = "blank")
for(i in 1:3){
	for(j in 1:3){
		if(i>j){
				aa <- putPlot(aa,ggplot(data.frame(x=dt.pc.load[,j],y=dt.pc.load[,i]))+
						geom_point(aes(x=x,y=y)) + geom_density2d(aes(x=x,y=y)),i,j)
		}
	}
}
pdf("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/02.cluster//enrich.gsea.pv.tb.logp.pca_pairPlot.pdf")
print(aa)
dev.off()

library(rococo)
gk.cor <- function(x){
	gk.dt <- matrix(0,nr=NROW(x),nc=NROW(x))
	dimnames(gk.dt) <- list(rownames(x),rownames(x))
	for(i in 1:NROW(x)){
		for(j in 1:NROW(x)){
			gk.dt[i,j] <- rococo(as.numeric(x[i,]),as.numeric(x[j,]),similarity="classical")
		}
	}
	return(gk.dt)
}

#abbreviations
#distance: euclidean (EUC), manhattan (MAN), maximum (SUP), canberra (CAN), binary (BIN), minkowski (MIN)
#cor: pearson (PE), kendall (KE), spearman (SP), Goodman-Kruskal (GK), cosine (COS), 
#cluster: kmean (KMN), k-medoid (KM), hcl_average (AL), hcl_complete (CL), hcl_single(SL), hcl_ward.D (WDL),
#hcl_ward.D2 (WD2L), hcl_mcquitty (MCL), hcl_median (MEL), hcl_centroid (CEL)
abr <- c("euclidean","manhattan","maximum","canberra","binary","minkowski","pearson","kendall","spearman","Goodman-Kruskal","cosine")
names(abr) <- c("EUC","MAN","SUP","CAN","BIN","MIN","PE","KE","SP","GK","COS")
abr <- c(abr,"kmean","k-medoid","hcl_average","hcl_complete","hcl_single","hcl_ward.D","hcl_ward.D2","hcl_mcquitty","hcl_median","hcl_centroid")
names(abr) <- c(names(abr)[1:11],"KMN","KM","AL","CL","SL","WDL","WD2L","MCL","MEL","CEL")
abr <- c(abr,"TOM");names(abr)[22] <- c("TOM")


#distance: euclidean (EUC), manhattan (MAN), maximum (SUP), canberra (CAN), binary (BIN), minkowski (MIN)
#cor: pearson (PE), kendall (KE), spearman (SP), Goodman-Kruskal (GK), cosine (COS), 

library(lsa)
#calculate distance
dst.lst <- list()
dst.lst[["EUC"]] <- dist(dt,method="euclidean")
dst.lst[["MAN"]] <- dist(dt,method="manhattan")
dst.lst[["SUP"]] <- dist(dt,method="maximum")
dst.lst[["CAN"]] <- dist(dt,method="canberra")
dst.lst[["BIN"]] <- dist(dt,method="binary")
dst.lst[["MIN"]] <- dist(dt,method="minkowski")
dst.lst[["PE"]] <- as.dist(1-cor(t(dt),method="pearson"))
dst.lst[["KE"]] <- as.dist(1-cor(t(dt),method="kendall"))
dst.lst[["SP"]] <- as.dist(1-cor(t(dt),method="spearman"))
dst.lst[["GK"]] <- as.dist(1-gk.cor(dt))
dst.lst[["COS"]] <- as.dist(1-cosine(t(dt)))
if(identical("y","y")){
	load("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/01.wgcna/enrich.gsea.pv.tb.logp.networkConstruction.TOM.RData")
	dst.lst[["TOM"]] <- as.dist(dissTOM)
}
save(dst.lst,file="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/02.cluster//enrich.gsea.pv.tb.logp.dst.RData")

#plot multidimentional scaling
mds.dt <- data.frame()
for(i in names(dst.lst)){
	tmp.dt <- cmdscale(dst.lst[[i]])
	mds.dt <- rbind(mds.dt,tmp.dt)
}
mds.dt <- cbind(mds.dt,rep(abr[names(dst.lst)],each=nrow(dt)))
colnames(mds.dt) <- c("ScalingDimension1","ScalingDimension2","Method")
pdf("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/02.cluster//enrich.gsea.pv.tb.logp.dst.mds_plot.pdf")
aa <- ggplot(mds.dt) + geom_point(aes(x=ScalingDimension1,y=ScalingDimension2))+
		geom_density2d(aes(x=ScalingDimension1,y=ScalingDimension2))+
		facet_wrap(~Method,scales="free")
print(aa)		
dev.off()


#cluster: kmean (KMN), k-medoid (KM), hcl_average (AL), hcl_complete (CL), hcl_single(SL), hcl_ward.D (WDL),
#hcl_ward.D2 (WD2L), hcl_mcquitty (MCL), hcl_median (MEL), hcl_centroid (CEL)

#hierarchical clustering
hcl.lst <- list()
for(i in names(dst.lst)){
	hcl.lst[[i]][["WDL"]] <- hclust(dst.lst[[i]],method="ward.D")
	hcl.lst[[i]][["WD2L"]] <- hclust(dst.lst[[i]],method="ward.D2")
	hcl.lst[[i]][["AL"]] <- hclust(dst.lst[[i]],method="average")
	hcl.lst[[i]][["CL"]] <- hclust(dst.lst[[i]],method="complete")
	hcl.lst[[i]][["SL"]] <- hclust(dst.lst[[i]],method="single")
	hcl.lst[[i]][["MCL"]] <- hclust(dst.lst[[i]],method="mcquitty")
	hcl.lst[[i]][["MEL"]] <- hclust(dst.lst[[i]],method="median")
	hcl.lst[[i]][["CEL"]] <- hclust(dst.lst[[i]],method="centroid")
}

save(hcl.lst,dst.lst,abr,file="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/02.cluster//enrich.gsea.pv.tb.logp.hcl.RData")

#cophenetic correlation
cop.mt <- matrix(0,nr=length(hcl.lst[[1]]),nc=length(hcl.lst))
dimnames(cop.mt) <- list(x=abr[names(hcl.lst[[1]])],y=abr[names(hcl.lst)])
for(i in names(dst.lst)){
	for(j in names(hcl.lst[[1]])){
		cop.mt[abr[j],abr[i]] <- cor(dst.lst[[i]],cophenetic(hcl.lst[[i]][[j]]))
	}
}
write.csv(cop.mt,file="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/02.cluster//enrich.gsea.pv.tb.logp.cophenetic_cor.csv")

#dynamic cut tree
#here we using dynamicTreeCut to select the cluster number, and then evaluate various 
#clustering methods
library(dynamicTreeCut)
library(clusterSim)
dynamic.tree <- list()
hcl.lst.al <- list()
dynamic.tree.idx <- matrix(nc=length(hcl.lst),nr=5)
dimnames(dynamic.tree.idx) <- list(c("DB","CH","Gamma","CI","S"),abr[names(hcl.lst)])
for(i in names(hcl.lst)){
	dynamic.tree[[i]] <- cutreeDynamic(dendro=hcl.lst[[i]][["AL"]],minClusterSize=2,distM=as.matrix(dst.lst[[i]]),
			                       cutHeight=as.numeric(quantile(hcl.lst[[1]][["AL"]]$height,0.90)))
	names(dynamic.tree[[i]]) <- hcl.lst[[i]][["AL"]]$labels
	if(length(which(dynamic.tree[[i]]==0)) > length(dynamic.tree[[i]])*0.5){
		dynamic.tree.idx["DB",abr[i]] <- "NA"
		dynamic.tree.idx["CH",abr[i]] <- "NA"
		dynamic.tree.idx["Gamma",abr[i]] <- "NA"
		dynamic.tree.idx["CI",abr[i]] <- "NA"
		dynamic.tree.idx["S",abr[i]] <- "NA"
		dynamic.tree <- dynamic.tree[-pmatch(i,names(dynamic.tree))]
	}else{
		dynamic.tree.idx["DB",abr[i]] <- index.DB(x=dt,cl=dynamic.tree[[i]],centrotypes="centroids")$DB #Minimum value of the index
		dynamic.tree.idx["CH",abr[i]] <- index.G1(x=as.data.frame(dt),cl=dynamic.tree[[i]]) #Maximum value of the index
		dynamic.tree.idx["Gamma",abr[i]] <- index.G2(d=dst.lst[[i]],cl=dynamic.tree[[i]]) #Maximum value of the index
		dynamic.tree.idx["CI",abr[i]] <- index.G3(d=dst.lst[[i]],cl=dynamic.tree[[i]]) #Minimum value of the index
		dynamic.tree.idx["S",abr[i]] <- index.S(d=dst.lst[[i]],cl=dynamic.tree[[i]]) #Maximum value of the index
		hcl.lst.al[[i]] <- hcl.lst[[i]][["AL"]]
	}
}
save(dynamic.tree,hcl.lst.al,abr,dst.lst,file="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/02.cluster//enrich.gsea.pv.tb.logp.tree.RData")

#select optimized methods
dynamic.tree.idx.slc <- vector("character",length=dim(dynamic.tree.idx)[1])
names(dynamic.tree.idx.slc) <- rownames(dynamic.tree.idx)
dynamic.tree.idx.slc["DB"] <- colnames(dynamic.tree.idx)[which.min(dynamic.tree.idx["DB",])]
dynamic.tree.idx.slc["CH"] <- colnames(dynamic.tree.idx)[which.max(dynamic.tree.idx["CH",])]
dynamic.tree.idx.slc["Gamma"] <- colnames(dynamic.tree.idx)[which.max(dynamic.tree.idx["Gamma",])]
dynamic.tree.idx.slc["CI"] <- colnames(dynamic.tree.idx)[which.min(dynamic.tree.idx["CI",])]
dynamic.tree.idx.slc["S"] <- colnames(dynamic.tree.idx)[which.max(dynamic.tree.idx["S",])]

write.csv(data.frame(dynamic.tree.idx,Optimize=dynamic.tree.idx.slc),file="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/02.cluster//enrich.gsea.pv.tb.logp.tree_idx.csv")


#cluster: kmean (KMN), k-medoid (KM), K-centroids (KC)
#select cluster numbers using gap statistics
library(cluster)
kn <- NROW(dt)/2
gap.kn <- clusGap(dt,FUN=pam,K.max=kn,B=100)
pdf("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/02.cluster//enrich.gsea.pv.tb.logp.gap_stat.pdf")
plot(gap.kn)
dev.off()
sink("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/02.cluster//enrich.gsea.pv.tb.logp.gap_stat.txt")
print(gap.kn,method="globalmax")
print(gap.kn,method="firstmax")
print(gap.kn,method="Tibs2001SEmax")
print(gap.kn,method="firstSEmax")
print(gap.kn,method="globalSEmax")
sink()

#after determining the cluster number, then fuzzy clustering can be used
#fanny {cluster}

for(i in names(dst.lst)){
	clst.lst[[i]][["KMN"]] <- kmean(dt,centers=kn,nstart=25)
	clst.lst[[i]][["KM"]] <- pam(dt,k=kn) #from cluster package
	clst.lst[[i]][["KC"]] <- kcca #from flexclust package
}

#clusterRepro
#pvclust

#do model based clustering
#et.mc <- Mclust(et.sig.log,G=1:20)
#plot(et.mc,what = "BIC")


q('no')
