
date()
library(piano)
#read in the gene level statistics
tmp = read.table("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/02.Enrich/00.gstat//G5_VS_S5",stringsAsFactors=FALSE)
glevel = tmp[,2]
names(glevel) = tmp[,1]
direct = tmp[,4]
names(direct) = tmp[,1]

#read in gene sets data
tmp = scan("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/01.prep/Kappa/gSets.gmt.filt.gmt.filtkappa.gmt",what="character",sep="\n")
tmp2 = sapply(tmp,function(x){
		x1 = strsplit(x,"\t")[[1]]
		x2 = paste(x1[-c(1,2)],x1[1],sep="\t")
		x2
		})
tmp2 = unlist(tmp2)
gene2gensets = data.frame(gene=sapply(tmp2,function(x) strsplit(x,"\t")[[1]][1]),
		set=sapply(tmp2,function(x) strsplit(x,"\t")[[1]][2]))
row.names(gene2gensets) <- 1:length(tmp2)
gscs = loadGSC(gene2gensets)
if(length(grep("fisher","gsea"))==1){
	gsaRes.fisher = runGSA(geneLevelStats=glevel,directions=direct,geneSetStat="fisher",
			signifMethod="geneSampling",adjMethod="fdr",gsc=gscs,gsSizeLim=c(4,3000),ncpus=4)
		GSAsummaryTable(gsaRes.fisher,save=TRUE,file="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/02.Enrich//GSEA//G5_VS_S5.enrich/G5_VS_S5.piano.fisher") 
}

if(length(grep("stouffer","gsea"))==1){
	gsaRes.stouffer = runGSA(geneLevelStats=glevel,directions=direct,geneSetStat="stouffer",
			signifMethod="geneSampling",adjMethod="fdr",gsc=gscs,gsSizeLim=c(4,3000),ncpus=4)
		GSAsummaryTable(gsaRes.stouffer,save=TRUE,file="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/02.Enrich//GSEA//G5_VS_S5.enrich/G5_VS_S5.piano.stouffer")
}
 
if(length(grep("reporter","gsea"))==1){
	gsaRes.reporter = runGSA(geneLevelStats=glevel,directions=direct,geneSetStat="reporter",
			signifMethod="geneSampling",adjMethod="fdr",gsc=gscs,gsSizeLim=c(4,3000),ncpus=4)
		GSAsummaryTable(gsaRes.reporter,save=TRUE,file="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/02.Enrich//GSEA//G5_VS_S5.enrich/G5_VS_S5.piano.reporter")
}

if(length(grep("tailStrength","gsea"))==1){
	gsaRes.tailStrength = runGSA(geneLevelStats=glevel,directions=direct,geneSetStat="tailStrength",
			signifMethod="geneSampling",adjMethod="fdr",gsc=gscs,gsSizeLim=c(4,3000),ncpus=4)
		GSAsummaryTable(gsaRes.tailStrength,save=TRUE,file="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/02.Enrich//GSEA//G5_VS_S5.enrich/G5_VS_S5.piano.tailStrength")
}

if(length(grep("wilcoxon","gsea"))==1){
	gsaRes.wilcoxon = runGSA(geneLevelStats=glevel,directions=direct,geneSetStat="wilcoxon",
			signifMethod="geneSampling",adjMethod="fdr",gsc=gscs,gsSizeLim=c(4,3000),ncpus=4)
		GSAsummaryTable(gsaRes.wilcoxon,save=TRUE,file="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/02.Enrich//GSEA//G5_VS_S5.enrich/G5_VS_S5.piano.wilcoxon")
}

if(length(grep("page","gsea"))==1){
	gsaRes.page = runGSA(geneLevelStats=direct,geneSetStat="page",
			signifMethod="geneSampling",adjMethod="fdr",gsc=gscs,gsSizeLim=c(4,3000),ncpus=4)
		GSAsummaryTable(gsaRes.page,save=TRUE,file="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/02.Enrich//GSEA//G5_VS_S5.enrich/G5_VS_S5.piano.page")
}

if(length(grep("gsea","gsea"))==1){
	gsaRes.gsea = runGSA(geneLevelStats=direct,geneSetStat="gsea",
			signifMethod="geneSampling",adjMethod="fdr",gsc=gscs,gsSizeLim=c(4,3000),ncpus=4)
		GSAsummaryTable(gsaRes.gsea,save=TRUE,file="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/02.Enrich//GSEA//G5_VS_S5.enrich/G5_VS_S5.piano.gsea")
}

if(length(grep("maxmean","gsea"))==1){
	gsaRes.maxmean = runGSA(geneLevelStats=direct,geneSetStat="maxmean",
			signifMethod="geneSampling",adjMethod="fdr",gsc=gscs,gsSizeLim=c(4,3000),ncpus=4)
		GSAsummaryTable(gsaRes.maxmean,save=TRUE,file="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/02.Enrich//GSEA//G5_VS_S5.enrich/G5_VS_S5.piano.maxmean")
}

#save(gsaRes.fisher,gsaRes.gsea,gsaRes.maxmean,gsaRes.page,gsaRes.reporter,gsaRes.stouffer,gsaRes.tailStrength,file="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/02.Enrich//GSEA//G5_VS_S5.enrich/G5_VS_S5.RData")

#consensus
#resList <- list(gsaRes.fisher,gsaRes.gsea,gsaRes.maxmean,gsaRes.page,gsaRes.reporter,gsaRes.stouffer,gsaRes.tailStrength)
#names(resList) <- c("fisher","gsea","maxmean","page","reporter","stouffer","tailStrength")
#ch <- consensusHeatmap(resList,cutoff=30,method="mean")
date()
q('no')
