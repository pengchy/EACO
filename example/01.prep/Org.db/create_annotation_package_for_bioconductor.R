
library(AnnotationForge)
options(stringsAsFactors=FALSE)
fChr <- read.table("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/00.data/geneid2scfid",sep="\t")[,1:2]
colnames(fChr) <- c("GID","CHROMOSOME")

fGO <- read.table("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/00.data/MetaData.GO.gmt.3cl",sep="\t",quote = "")
fGO <- data.frame(fGO[,1:2],rep("IEA",dim(fGO)[1]))
colnames(fGO) <- c("GID","GO","EVIDENCE")

fSym <- read.table("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/00.data/gene.id.desc",sep="\t",quote = "")
fSym <- data.frame(fSym[,1],fSym[,1],fSym[,2])
fSym <- fSym[which(!is.na(fSym[,3])),]
colnames(fSym) <- c("GID","SYMBOL","GENENAME")

makeOrgPackage(gene_info=fSym, chromosome=fChr, go=fGO,
		version="0.1",
		maintainer="Pengcheng Yang <yangpc@mail.biols.ac.cn>",
		author="Pengcheng Yang <yangpc@mail.biols.ac.cn>",
		outputDir = "/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/01.prep/Org.db",
		tax_id="7004",
		genus="Locusta",
		species="migratoria",
		goTable="go")

packge.name <- paste("org.",substr("Locusta",1,1),"migratoria",".eg.db",sep="")
pckg <- paste("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/01.prep/Org.db/",package.name,sep="")
#install the packages
install.packages(pckg,repos=NULL)

makeDBPackage("GO_DB",
		affy=FALSE,
		prefix="lmi_GO",
		fileName="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/00.data/MetaData.GO.gmt.3cl",
		baseMapType="gb",
		outputDir = tmpout,
		version="1.0.0",
		manufacturer = "Affymetrix",
		chipName = "Human Cancer G110 Array",
		manufacturerUrl = "http://www.affymetrix.com")

hcg110_IDs = system.file("extdata",
		"hcg110_ID",
		package="AnnotationDbi")
tmpout = "./"
	makeDBPackage("HUMANCHIP_DB",
			affy=FALSE,
			prefix="hcg110",
			fileName=hcg110_IDs,
			baseMapType="gb",
			outputDir = tmpout,
			version="1.0.0",
			manufacturer = "Affymetrix",
			chipName = "Human Cancer G110 Array",
			manufacturerUrl = "http://www.affymetrix.com")

GO_DB = createSimpleBimap(
		"go",
		"GID",
		"GO",
		org.Lmigratoria.eg.db:::datacache,
		"GO_DB",
		"org.Lmigratoria.eg.db"
		)
