
date()
library(GOstats) 
library("GSEABase")
options(stringsAsFactors=FALSE)
tmp <- scan("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/00.data/MetaData.GO.gmt",what="character",sep="\n")
gsets <- sapply(tmp,function(x) strsplit(x,"\t")[[1]][-c(1,2)])
names(gsets) <- sapply(tmp,function(x) strsplit(x,"\t")[[1]][1])
setDes <- sapply(tmp,function(x) strsplit(x,"\t")[[1]][2])
names(setDes) <- names(gsets)

if(identical("N","N")){
	universe <- as.character(unlist(gsets))
}else{
	universe <- scan("N",sep="\n",what="character")
}

suppid <- scan("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/02.Enrich/00.gstat/Gadult-VS-Sadult.down",sep="\n",what="character")

if(identical("GO","GO")){
	goframeData <- data.frame(frame.go_id=rep(names(gsets),sapply(gsets,length)),
			frame.Evidence=rep("IEA",length(unlist(gsets))),
			frame.gene_id=as.character(unlist(gsets)))
	goFrame=GOFrame(goframeData,organism="Locusta migratoria")
	goAllFrame=GOAllFrame(goFrame)

	gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
	params <- GSEAGOHyperGParams(name="Custom GSEA based annot Params",
			geneSetCollection=gsc,
			geneIds = suppid,
			universeGeneIds = universe,
			ontology = "MF",
			pvalueCutoff = 0.05,
			conditional = FALSE,
			testDirection = "over")
	Over <- hyperGTest(params)

}


