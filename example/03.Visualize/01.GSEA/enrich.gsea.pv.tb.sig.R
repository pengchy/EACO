
date()
#filtering the enrichment results
options(stringsAsFactors=FALSE)
et <- read.table("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/02.Enrich/GSEA/enrich.gsea.pv.tb",header=TRUE,sep="\t",check.names=FALSE)
rownames(et) <- toupper(et[,1])
et <- as.matrix(et[,-1])
et[is.na(et)] <- 1
et[et==0] <- 1e-5
etcfg <- read.table("/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/02.Enrich/GSEA/Etb.cfg",sep="\t",check.names=FALSE)
et <- et[,etcfg[,1]]
#filt enrichment result
et.sig <- et[apply(et,1,function(x) any(x<0.05)),]
et.sig.log <- -log10(et.sig)
write.table(et.sig.log,file="/panfs/home/kang/yangpc/bin/EnrichPipeline/EACO_r20150201/example/03.Visualize/01.GSEA/enrich.gsea.pv.tb.logp",sep="\t",quote=FALSE)

date()
q('no')
