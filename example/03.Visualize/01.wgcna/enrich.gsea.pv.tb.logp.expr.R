
#reading data
library(WGCNA) 
options(stringsAsFactors = FALSE)
expr <- read.table("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/enrich.gsea.pv.tb.logp",check.names=FALSE)
expr <- expr[apply(expr,1,function(x) any(x>0)),]
#output the filtered expression value
write.table(expr,file="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/01.wgcna/enrich.gsea.pv.tb.logp.after_filt",sep="\t",quote=FALSE)
datExpr0 <- as.data.frame(t(expr))
if(identical("y","n")){
	gsg <- goodSamplesGenes(expr,verbos=3)
	if (!gsg$allOK){
		if (sum(!gsg$goodGenes)>0)
			printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
		if (sum(!gsg$goodSamples)>0)
			printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
		datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
	}
}

datExpr0.dimn <- dimnames(datExpr0)
datExpr0 <- log2(datExpr0+1)
datExpr0 <- t(scale(t(datExpr0),center=TRUE,scale=TRUE))
datExpr0 <- scale(datExpr0,center=TRUE,scale=TRUE)
dimnames(datExpr0) <- datExpr0.dimn
#checking data quality
if(identical("y","n")){
	sampleTree = flashClust(dist(datExpr0), method = "average");
	pdf("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/01.wgcna/enrich.gsea.pv.tb.logp.Sample_clustering_to_detect_outliers.pdf",height=9,width=12)
	par(cex = 0.6);
	par(mar = c(0,4,2,0))
	plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
		cex.axis = 1.5, cex.main = 2)
	dev.off()
}
datExpr <- datExpr0
#construct coexpression network
enableWGCNAThreads()
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
#select beta
pdf("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/01.wgcna/enrich.gsea.pv.tb.logp.Power_selection.pdf",height=5,width=9)
par(mfrow = c(1,2));
cex1 = 0.9; 
# Scale-free topology fit index as a function of the soft-thresholding power 
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
save(datExpr,expr,sft,file="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/01.wgcna/enrich.gsea.pv.tb.logp.expr.RData")
q('no')
