
library(WGCNA)
load("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/01.wgcna/enrich.gsea.pv.tb.logp.expr.RData")
#construct network step-by-step
enableWGCNAThreads()
#Co-expression similarity and adjacency
adjacency <- adjacency(datExpr, power = 1,type="signed");
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dimnames(TOM) <- dimnames(adjacency)
dissTOM = 1-TOM
save(adjacency,dissTOM,file="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/01.wgcna/enrich.gsea.pv.tb.logp.networkConstruction.TOM.RData")
geneTree = flashClust(as.dist(dissTOM), method = "average")
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
		deepSplit = 2, pamRespectsDendro = FALSE,
		minClusterSize = 20);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and the module colors underneath
pdf("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/01.wgcna/enrich.gsea.pv.tb.logp.DendroAndColors.pdf")
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
		dendroLabels = FALSE, hang = 0.03,
		addGuide = TRUE, guideHang = 0.05,
		main = "Gene dendrogram and module colors")
dev.off()

if(identical("N","Y")){
	merged <- mergeCloseModules(datExpr, dynamicColors, cutHeight = 0.2,verbos=3)
	pdf("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/01.wgcna/enrich.gsea.pv.tb.logp.DendroAndColors.after_merge.pdf")
	plotDendroAndColors(geneTree, cbind(dynamicColors,merged$colors),
			c("Dynamic Tree Cut","Merged Dynamic"),dendroLabels = FALSE, hang = 0.03,
			addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
	dev.off()
	dynamicColors <- merged$colors
}

save(datExpr,dynamicColors,adjacency,TOM,geneTree,file="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/01.wgcna/enrich.gsea.pv.tb.logp.for_net_trait.RData")

save(datExpr,expr,dynamicMods,dynamicColors,geneTree,
			file = "/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150101/example/03.Visualize/01.wgcna/enrich.gsea.pv.tb.logp.networkConstruction.RData")
q('no')
